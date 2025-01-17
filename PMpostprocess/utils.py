import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import fftconvolve
from scipy.optimize import minimize
import pickle
from pprint import pprint


def wire_convolution_kernel(x, r):
    """
    Generate a convolution kernel for a cylindrical wire of radius r.
    The kernel assumes the wire cross-section is circular and projects it onto the beam plane.
    """
    x = x-np.mean(x)
    kernel = np.where(np.abs(x) <= r, np.sqrt(r**2 - x**2), 0)  # Circle equation
    kernel /= np.sum(kernel)  # Normalize kernel to ensure energy conservation
    return kernel


def gaussian_convolution_kernel(x, sigma):
    x = x-np.mean(x)
    kernel = np.exp(- 0.5*x**2/sigma**2)  # Circle equation
    kernel /= np.sum(kernel)  # Normalize kernel to ensure energy conservation
    return kernel


def convolve(signal, kernel):
    full_convolved = fftconvolve(signal, kernel, mode='full')
    # Trim to match the original size
    original_size = len(signal)
    kernel_size = len(kernel)
    start = (kernel_size - 1) // 2
    end = start + original_size
    return 0.5*(full_convolved[start+1:end+1]+full_convolved[start:end])
    
    
def exponential_decay(iteration, initial_strength=0.9, decay_rate=0.1):
    return initial_strength * np.exp(-decay_rate * iteration)
    
    
def richardson_lucy_deconvolution(signal, kernel, iterations=20, reg_kernel=None, decay_function=None):
    """
    Perform Richardson-Lucy deconvolution.
    Args:
        signal: Measured signal (1D array).
        kernel: Convolution kernel (1D array).
        iterations: Number of iterations for deconvolution.
    Returns:
        Deconvolved signal (1D array).
    """
    minval = np.min(signal)
    maxval = np.max(signal)
    eps = 1e-3
    s = (signal + minval)/(maxval-minval) + eps
    kernel_flipped = kernel[::-1]  # Flip kernel for convolution
    estimated = 0.5*np.ones_like(s)

    for i in range(iterations):
        convolved = convolve(estimated, kernel)
        ratio = (s + eps) / (convolved + eps)  # Prevent division by zero
        estimated *= convolve(ratio, kernel_flipped) #fftconvolve(ratio, kernel_flipped, mode='same')  # Update estimate
        if reg_kernel is not None:
            decay_factor = decay_function(i) if decay_function is not None else exponential_decay(i)
            smoothed = convolve(estimated, reg_kernel)
            estimated = (1 - decay_factor) * estimated + decay_factor * smoothed
    
    return (estimated-eps)*(maxval-minval) - minval
    
def compute_weights(positions):
    """
    Compute weights for data fidelity and TV regularization based on positions.
    """
    spacing = np.diff(positions)
    weights = np.zeros(len(positions))
    reg_weights = np.zeros(len(positions) - 1)  # Regularization weights are between points

    # Data fidelity weights (average spacing for interior points)
    weights[1:-1] = 0.5 * (spacing[:-1] + spacing[1:])
    weights[0] = spacing[0]
    weights[-1] = spacing[-1]
    weights /= np.sum(weights)  # Normalize

    # Regularization weights (spacing directly between points)
    reg_weights[:] = spacing
    reg_weights /= np.sum(reg_weights)  # Normalize

    return weights, reg_weights

def tv_l1_loss_weighted(estimated, signal, convolved, weights, reg_weights, lambda_reg):
    """
    Compute the loss function with weighted TV regularization and weighted L1 data fidelity.
    """
    # Weighted L1 data fidelity term
    data_fidelity = np.sum(weights * np.abs(signal - convolved))
    
    # Weighted TV regularization term
    tv_regularization = lambda_reg * np.sum(np.abs(np.diff(estimated)) * reg_weights)
    
    return data_fidelity + tv_regularization

def gradient_tv_l1_loss_weighted(estimated, signal, convolved, kernel, weights, reg_weights, lambda_reg):
    """
    Compute the gradient of the loss function with weighted TV regularization
    and weighted L1 data fidelity.
    """
    # Weighted L1 data fidelity gradient
    fidelity_grad = -convolve(weights * np.sign(signal - convolved), kernel[::-1])
    
    # Weighted TV gradient (vectorized)
    tv_grad = np.zeros_like(estimated)
    tv_grad[:-1] += lambda_reg * np.sign(estimated[:-1] - estimated[1:]) * reg_weights
    tv_grad[ 1:] -= lambda_reg * np.sign(estimated[:-1] - estimated[1:]) * reg_weights
    
    return fidelity_grad + tv_grad


def loss_and_grad(x, signal, kernel, weights, reg_weights, lambda_reg):
    """
    Combined loss and gradient computation for the optimizer.
    """
    # Compute convolution once
    convolved = convolve(x, kernel)

    # Compute loss
    loss = tv_l1_loss_weighted(x, signal, convolved, weights, reg_weights, lambda_reg)

    # Compute gradient
    grad = gradient_tv_l1_loss_weighted(x, signal, convolved, kernel, weights, reg_weights, lambda_reg)

    return loss, grad


def deconvolve_signal(positions, signal, kernel, initial_guess=None, initial_lambda_reg=0.1, decay_rate=0.01, max_iters=100):
    """
    Deconvolve the given signal using TV-L1 regularization with decaying lambda.
    
    Args:
        positions: Array of positions corresponding to the signal.
        signal: Measured signal (1D array).
        kernel: Convolution kernel (1D array).
        initial_lambda_reg: Initial regularization strength for TV.
        decay_rate: Rate at which lambda_reg decays per iteration.
        max_iters: Maximum number of iterations.
        
    Returns:
        Deconvolved signal (1D array).
    """
    # Compute weights based on positions
    weights, reg_weights = compute_weights(positions)
    
    # Use an initial guess based on the median of the signal
    if initial_guess is None:
        initial_guess = np.full_like(signal, np.median(signal))
    
    # Track lambda_reg and iteration count
    lambda_reg = initial_lambda_reg
    iteration = 0

    def callback(xk):
        """
        Callback function to decay lambda_reg over iterations.
        """
        nonlocal lambda_reg, iteration
        iteration += 1
        lambda_reg = max(0, initial_lambda_reg * (1 - decay_rate * iteration))

    # Optimization using SciPy minimize
    result = minimize(
        fun=lambda x: loss_and_grad(x, signal, kernel, weights, reg_weights, lambda_reg),
        x0=initial_guess,
        jac=True,  # Indicate that `fun` also provides the gradient
        method='L-BFGS-B',
        bounds=[(0, None)] * len(signal),  # Enforce positivity
        callback=callback,
        options={'maxiter': max_iters}
    )

    # Return the optimized signal
    return result
    
def get_istart_iend_profile(y, threshold=None):
    argmax = np.argmax(y)
    if threshold is None:
        threshold = 0.02 * y[argmax]  # 5% of the maximum value
    
    # Find start index using NumPy searchsorted (efficient linear search)
    start_candidates = np.where(y[:argmax] <= threshold)[0]
    start_index = start_candidates[-1] if start_candidates.size > 0 else 0
    while start_index > 0 and 0<y[start_index] < y[start_index + 1]:
        start_index -= 1
    
    # Find end index
    end_candidates = np.where(y[argmax:] <= threshold)[0] + argmax
    end_index = end_candidates[0] if end_candidates.size > 0 else len(y) - 1
    while end_index < len(y) - 1 and 0<y[end_index] < y[end_index - 1]:
        end_index += 1
    
    return start_index, end_index

def estimate_noise_floor(start_index, end_index, y):
    noise_region = np.concatenate([y[:start_index], y[end_index + 1:]])
    if len(noise_region) > 10:
        # Directly calculate mean and std without additional steps
        noise_offset = np.mean(noise_region)
        noise_std = np.std(noise_region)
    else:
        noise_offset = 0
        noise_std = None    
    return noise_offset, noise_std    

def cut_boundary(y, start_index, end_index, offset):
    y = np.maximum(y - offset, 0)
    y[:start_index] = 0
    y[end_index:] = 0
    return y
    
def gaussian_smooth_with_deconvolution(x, y, scale_factor=2.0, edge_cut=True):
    spacing = np.diff(x)
    # Estimate global sigma
    sigma = scale_factor * np.mean(spacing)
    kernel = gaussian_convolution_kernel(x, sigma)
    reg_kernel = gaussian_convolution_kernel(x, 0.3*sigma)
    # deconv_kernel = gaussian_convolution_kernel(x, (1+0.04)**0.5*sigma)
    deconv_kernel = kernel
    convolved = convolve(y, kernel)  # Forward convolution
    start_index, end_index = get_istart_iend_profile(convolved)
    if edge_cut:
        #convolved = cut_boundary(convolved, start_index, end_index, 0)
        convolved = cut_boundary(convolved,start_index,end_index,0)
    deconvolved = richardson_lucy_deconvolution(convolved, deconv_kernel, reg_kernel = reg_kernel)
    if finetune_deconv:
        result = deconvolve_signal(x, convolved, kernel, initial_guess=deconvolved, initial_lambda_reg=0.1)
        deconvolved = result.x
    if edge_cut:
        start_index, end_index = get_istart_iend_profile(deconvolved)
        deconvolved = cut_boundary(deconvolved,start_index,end_index,0)
    
    return convolved, deconvolved, start_index, end_index
    
def measure_rms_size(positions,signal):
    positions_mean = np.sum(positions * signal) / np.sum(signal)
    shifted_positions = positions - positions_mean
    return np.sqrt(np.sum(signal * shifted_positions**2) / np.sum(signal))
    
def kde_estimate(signal, positions, bandwidth=1.0):
    # Apply Kernel Density Estimation
    kde = gaussian_kde(positions, weights=signal, bw_method=bandwidth / np.std(positions))
    kde_values = kde(positions)  # Evaluate the KDE at the positions    
    return kde_values

def scale_adjustment_factor(kde_est, processed_signal):
    numerator = np.sum(kde_est * processed_signal)
    denominator = np.sum(kde_est ** 2)
    return numerator / denominator