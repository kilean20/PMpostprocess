from collections.abc import Iterable
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import fftconvolve
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import pickle
from pprint import pprint
import time
import logging

logger = logging.getLogger(__name__)

def wire_convolution_kernel(x, r):
    """
    Generate a convolution kernel for a cylindrical wire of radius r.
    The kernel assumes the wire cross-section is circular and projects it onto the beam plane.
    """
    x = x - np.mean(x)  # Center x around its mean
    # Ensure no negative values are passed to sqrt by clamping
    kernel = np.where(np.abs(x) <= r, np.sqrt(np.maximum(0, r**2 - x**2)), 0)  # Circle equation
    kernel /= np.sum(kernel)  # Normalize kernel to ensure energy conservation
    return kernel



def gaussian_convolution_kernel(x, sigma):
    if sigma < 1e-10:
        kernel = np.ones_like(x)
    else:
        x = x-np.mean(x)
        kernel = np.exp(- 0.5*x**2/sigma**2)  # Circle equation
    kernel /= np.sum(kernel)  # Normalize kernel to ensure energy conservation
    return kernel


def convolve(y, kernel):
    full_convolved = fftconvolve(y, kernel, mode='full')
    # Trim to match the original size
    original_size = len(y)
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
            smoothed = convolve(estimated, reg_kernel)
            decay_factor = decay_function(i) if decay_function is not None else exponential_decay(i)
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


def deconvolve_fit(positions, signal, kernel, initial_guess=None, initial_lambda_reg=0.1, decay_rate=0.01, max_iters=100):
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
    
def measure_rms_size(positions,signal):
    positions_mean = np.sum(positions * signal) / np.sum(signal)
    shifted_positions = positions - positions_mean
    return np.sqrt(np.sum(signal * shifted_positions**2) / np.sum(signal))
    
    
def get_istart_iend_profile(y, threshold=None):
    y = y + np.quantile(y,0.02)
    argmax = np.argmax(y)
    if threshold is None:
        threshold = 0.02*y[argmax]  # 2% of the maximum value
    
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

def estimate_noise_floor(y,start_index, end_index):
    noise_floor_region = np.concatenate([y[:start_index], y[end_index + 1:]])
    if len(noise_floor_region) > 3:
        noise_offset = np.mean(noise_floor_region)
        noise_std = np.std(noise_floor_region)
    else:
        noise_offset = np.quantile(y,0.02)
        threshold = noise_offset + (noise_offset-np.min(y))
        noise_std = np.std(y[y<threshold])
    return noise_offset, noise_std

def cut_boundary(y, start_index, end_index, offset):
    y = np.maximum(y - offset, 0)
    y[:start_index] = 0
    y[end_index:] = 0
    y[y < 0] = 0
    return y
    
def next_power_of_2(n):
    if n < 2:
        return 2
    return 1 << n.bit_length()
        
def wire_deconvolution(x, y, r, wire_kernel_func=wire_convolution_kernel, is_x_uniform=False, extrapolate_factor=8, finetune_deconv=False):
    if is_x_uniform:
        xu = x
        yu = y
        n  = len(x)
    else:
        n = next_power_of_2(len(x))*extrapolate_factor
        xu = np.linspace(np.min(x),np.max(x),n)
        yu = interp1d(x, y)(xu)
    if r < np.diff(xu).mean():
        return xu, y
    wire_kernel = wire_kernel_func(xu, r)        
    reg_kernel = gaussian_convolution_kernel(xu, 0.3*r)
    deconvolved = richardson_lucy_deconvolution(yu, wire_kernel, reg_kernel = reg_kernel)
    if finetune_deconv:
        result = deconvolve_fit(xu, yu, wire_kernel, initial_guess=yu, initial_lambda_reg=0.1)
        deconvolved = result.x
    return xu, deconvolved

    
def gaussian_smooth_with_deconvolution(x, y, is_x_uniform=False, extrapolate_factor=8, edge_cut=True, finetune_deconv=False):

    spacing = np.diff(x)
    sigma = 2*np.mean(spacing)
    if is_x_uniform:
        xu = x
        yu = y
        n  = len(x)
    else:
        n = next_power_of_2(len(x))*extrapolate_factor
        xu = np.linspace(np.min(x),np.max(x),n)
        yu = interp1d(x, y)(xu)
    
    kernel = gaussian_convolution_kernel(xu, sigma)
    convolved = convolve(yu, kernel)
    istart, iend = get_istart_iend_profile(convolved)
    
    rms_tmp = measure_rms_size(xu[istart:iend],yu[istart:iend])
    if 0.15*rms_tmp < sigma:
        sigma = 0.15*rms_tmp
        kernel = gaussian_convolution_kernel(xu, sigma)
        convolved = convolve(yu, kernel)
        istart, iend = get_istart_iend_profile(convolved)
        
    reg_kernel = gaussian_convolution_kernel(xu, 0.3*sigma)
    deconv_kernel = kernel

    offset, std = estimate_noise_floor(yu,istart,iend)
    if edge_cut:
        convolved = cut_boundary(convolved,istart,iend,offset)
    deconvolved = richardson_lucy_deconvolution(convolved, deconv_kernel, reg_kernel = reg_kernel)
    if finetune_deconv:
        result = deconvolve_fit(xu, convolved, kernel, initial_guess=deconvolved, initial_lambda_reg=0.1)
        deconvolved = result.x
    if edge_cut:
        istart, iend = get_istart_iend_profile(deconvolved)
        deconvolved = cut_boundary(deconvolved,istart,iend,offset)
        
    #convolved = interp1d(xu, convolved)(x)
    #deconvolved = interp1d(xu, deconvolved)(x)
    #istart, iend = get_istart_iend_profile(convolved)
    #return convolved, deconvolved, istart, iend
    return xu, deconvolved
        
def estimate_two_noise_model(y,smoothed):
    scaler = np.max(smoothed)-np.min(smoothed)
    
    y = y/scaler
    smoothed = smoothed/scaler
    
    istart, iend = get_istart_iend_profile(smoothed)
    var_ysmooth = np.var(smoothed)
    var_noise = np.var(y-smoothed)
    
    noise_floor_offset, noise_floor_std = estimate_noise_floor(y,istart,iend)
    y = cut_boundary(y,istart,iend,noise_floor_offset)
    noise_signal_std = np.std(y[istart+1:iend-1]/smoothed[istart+1:iend-1])
    
    def loss(params):
        sigma1, sigma2 = params
        estimated_var_noise = sigma1**2 + sigma2**2 * var_ysmooth
        regularization = (sigma1/noise_floor_std-1)**2 + (sigma2/noise_signal_std-1)**2
        return 0.1*(estimated_var_noise/var_noise-1)**2 + 0.9*regularization

    result = minimize(loss, [noise_floor_std,noise_signal_std], bounds=[(0, 0.2), (0, 0.2)])
    noise_floor_std,noise_signal_std = result.x[0]*scaler,result.x[1]
    
    return noise_floor_std,noise_signal_std

def smooth_n_wire_deconvolve(x,y,r,extrapolate_factor=8):#,finetune_deconv=False):
    """
    x : position array
    y : signal array
    r : wire radius
    """    
    xu, smoothed = gaussian_smooth_with_deconvolution(x, y,is_x_uniform=False,extrapolate_factor=extrapolate_factor)
    yu = interp1d(x, y)(xu)
    n = len(xu)
    if r is not None and r > 0:
        xu, wire_deconvolved = wire_deconvolution(xu,smoothed,r,is_x_uniform=True)
    else:
        wire_deconvolved = smoothed
    
    istart, iend = get_istart_iend_profile(smoothed)
    noise_offset, _ = estimate_noise_floor(yu,istart,iend)
    yu = cut_boundary(yu,istart,iend,noise_offset)
    
    rms_beam_size   = measure_rms_size(xu, yu)
    rms_smooth_size = measure_rms_size(xu, smoothed)
    rms_deconv_size = measure_rms_size(xu, wire_deconvolved)
    
    smoothed = interp1d(xu,smoothed)(x)
    wire_deconvolved = interp1d(xu,wire_deconvolved)(x)
    
    return smoothed, wire_deconvolved, rms_beam_size, rms_smooth_size, rms_deconv_size
    
    
def rms_uncertainty_quantification(x,smoothed, noise_floor_std, noise_signal_std, r=None, nMC=32, extrapolate_factor = 2):
    arr_deconv_rms = np.zeros(nMC)
    arr_smooth_rms = np.zeros(nMC)
    arr_noisy_rms = np.zeros(nMC)
    arr_y_samples = np.zeros((nMC,len(x)))
    for i in range(nMC):
        y_ = smoothed*(1+noise_signal_std*np.random.randn(*x.shape)) + noise_floor_std*np.random.randn(*x.shape)
        _,_, noisy_rms, smooth_rms, deconv_rms = smooth_n_wire_deconvolve(x,y_,r,extrapolate_factor=extrapolate_factor)
        arr_smooth_rms[i] = smooth_rms
        arr_deconv_rms[i] = deconv_rms
        arr_noisy_rms[i] = noisy_rms
        arr_y_samples[i] = y_


    # arr_deconv_rms = arr_deconv_rms[np.logical_not(np.isnan(arr_deconv_rms))]
    # Calculating statistics
    # smooth_mean, smooth_std = np.nanmean(arr_smooth_rms), np.nanstd(arr_smooth_rms)
    # deconv_mean, deconv_std = np.nanmean(arr_deconv_rms), np.nanstd(arr_deconv_rms)
    # noisy_mean , noisy_std  = np.nanmean(arr_noisy_rms), np.nanstd(arr_noisy_rms)

    return {"arr_deconv_rms": arr_deconv_rms,
            "arr_y_samples": arr_y_samples,
            }

    
def process_profile_signal(x,y,r):
    smoothed, wire_deconvolved, rms_noisy, rms_smooth, rms_deconv = smooth_n_wire_deconvolve(x, y, r)
    noise_floor_std, noise_signal_std = estimate_two_noise_model(y, smoothed)
    stat = rms_uncertainty_quantification(x, smoothed, noise_floor_std, noise_signal_std, r, nMC=20)
    
    return {
        "smoothed": smoothed,
        "wire_deconvolved": wire_deconvolved,
        "rms_noisy": rms_noisy,
        "rms_smooth": rms_smooth,
        "rms_deconv": rms_deconv,
        "noise_floor_std": noise_floor_std,
        "noise_signal_std": noise_signal_std,
        "MC_stat": stat,
    }
    
    
def project_L3_to_xyuv(r6in, r12in1, r12in2, par_dict):
    r6in, r12in1, r12in2 = _convert_to_array(r6in, r12in1, r12in2)

    if not isinstance(r6in, np.ndarray):
        return _project_L3_to_xyuv_scalar(r6in, r12in1, r12in2, par_dict)
    else:
        return _project_L3_to_xyuv_vector(r6in, r12in1, r12in2, par_dict)

def _convert_to_array(*args):
    return [np.array(arg) if isinstance(arg, Iterable) and not isinstance(arg, np.ndarray) else arg for arg in args]

    

def _project_L3_to_xyuv_scalar(r6in, r12in1, r12in2, par_dict):
    coord = par_dict['coord']
    x = y = u = v = np.nan

    if coord.startswith("Luv"):
        u, v, primary = r6in, r12in1, r12in2 / math.sqrt(2)
        if coord == "Luvx":
            x, y = primary, _compute_sqrt(u**2 + v**2 - primary**2)
        elif coord == "Luvy":
            y, x = primary, _compute_sqrt(u**2 + v**2 - primary**2)

    elif coord.startswith("Lyx"):
        ang1, ang2 = par_dict['ang1'], par_dict['ang2']
        y, x = r6in * math.sin(math.radians(ang1)), r12in1 * math.cos(math.radians(ang2))
        if coord == "Lyxu":
            u, v = r12in2 * math.cos(math.radians(ang2 + 45)), _compute_sqrt(x**2 + y**2 - u**2)
        elif coord == "Lyxv":
            v, u = r12in2 * math.cos(math.radians(ang2 - 45)), _compute_sqrt(x**2 + y**2 - v**2)

    elif coord.startswith("S"):
        ang1 = par_dict['ang1']
        if coord == "Suxy":
            u, x, y = r6in, r12in1 * math.cos(math.radians(ang1)), r12in2 * math.sin(math.radians(ang1))
            v = _compute_sqrt(x**2 + y**2 - u**2)
        elif coord == "Sxyu":
            x, y, u = r6in * math.cos(math.radians(ang1)), r12in1 * math.sin(math.radians(ang1)), r12in2
            v = _compute_sqrt(x**2 + y**2 - u**2)
        elif coord == "Sbxy":
            x, y = r12in1 * math.cos(math.radians(ang1)), r12in2 * math.sin(math.radians(ang1))

    elif coord == "Fxy":
        x, y = r6in, r12in1


    cxy = np.nan   
    if x != 0 and y != 0:
        cxy = -(v**2 - u**2) / (2.0 * x * y)
    return x, y, cxy, u, v
    


def _project_L3_to_xyuv_vector(r6in, r12in1, r12in2, par_dict):
    assert len(r6in) == len(r12in1) == len(r12in2)
    coord = par_dict['coord']
    x, y, u, v = [np.full_like(r6in, np.nan) for _ in range(4)]  # Initialize all to np.nan

    if coord.startswith("Luv"):
        u, v, primary = r6in, r12in1, r12in2 / np.sqrt(2)
        if coord == "Luvx":
            x = primary
            mask = u**2 + v**2 - x**2 > 1e-10  # Add tolerance
            y[mask] = np.sqrt(u[mask]**2 + v[mask]**2 - x[mask]**2)
        elif coord == "Luvy":
            y = primary
            mask = u**2 + v**2 - y**2 > 1e-10
            x[mask] = np.sqrt(u[mask]**2 + v[mask]**2 - y[mask]**2)

    elif coord.startswith("Lyx"):
        ang1, ang2 = math.radians(par_dict['ang1']), math.radians(par_dict['ang2'])
        y, x = r6in * np.sin(ang1), r12in1 * np.cos(ang2)
        if coord == "Lyxu":
            u = r12in2 * np.cos(ang2 + math.pi / 4)
            v_mask = x**2 + y**2 - u**2 > 1e-10
            v[v_mask] = np.sqrt(x[v_mask]**2 + y[v_mask]**2 - u[v_mask]**2)
        elif coord == "Lyxv":
            v = r12in2 * np.cos(ang2 - math.pi / 4)
            u_mask = x**2 + y**2 - v**2 > 1e-10
            u[u_mask] = np.sqrt(x[u_mask]**2 + y[u_mask]**2 - v[u_mask]**2)

    elif coord.startswith("S"):
        ang1 = math.radians(par_dict['ang1'])
        if coord == "Suxy":
            x, y = r12in1 * np.cos(ang1), r12in2 * np.sin(ang1)
            u = r6in
            mask = x**2 + y**2 - u**2 > 1e-10
            v[mask] = np.sqrt(x[mask]**2 + y[mask]**2 - u[mask]**2)
        elif coord == "Sxyu":
            x, y = r6in * np.cos(ang1), r12in1 * np.sin(ang1)
            u = r12in2
            mask = x**2 + y**2 - u**2 > 1e-10
            v[mask] = np.sqrt(x[mask]**2 + y[mask]**2 - u[mask]**2)
        elif coord == "Sbxy":
            x, y = r12in1 * np.cos(ang1), r12in2 * np.sin(ang1)

    elif coord == "Fxy":
        x, y = r6in, r12in1

    # Compute cxy only for valid entries
    cxy = np.full_like(x, np.nan)
    valid_mask = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(u)) & (~np.isnan(v)) & (x != 0) & (y != 0)
    cxy[valid_mask] = -(v[valid_mask]**2 - u[valid_mask]**2) / (2.0 * x[valid_mask] * y[valid_mask])
    
    return x, y, cxy, u, v


def _compute_sqrt(value):
    return math.sqrt(value) if value > 0 else np.nan
    
    
#     fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=96)  # Create a figure with two subplots in a row
#     axes[0].plot(x, y, color="black", label=f"Noisy Profile,  $\\sigma_x$={rms_noisy:.3f} $\\pm$ {stat['noisy_std']:.3f} mm")
#     axes[0].plot(x, smoothed, color="green", label=f"Smoothed,   $\\sigma_x$={rms_smooth:.3f} $\\pm$ {stat['smooth_std']:.3f} mm")
#     axes[0].plot(x, wire_deconvolved, '--', color="red", label=f"DeConvolved, $\\sigma_x$={rms_deconv:.3f} $\\pm$ {stat['deconv_std']:.3f} mm")
#     axes[0].legend()
#     axes[0].set_xlabel("Position (mm)")
#     axes[0].set_ylabel("Signal Strength")
#     axes[0].set_ylim(0, 1.4 * max(y))
#     axes[0].set_title(f"Wire-thickness: {r:.3f} $mm$")

#     axes[1].plot(x, y - smoothed, 'k', label="Noisy - Smoothed")
#     for i in range(4):
#         y_ = stat['y_samples'][i]
#         axes[1].plot(x, y_ - smoothed, ':', label=f"Sample {i+1}")
#     axes[1].set_xlabel("Position (mm)")
#     axes[1].set_ylabel("Difference")
#     axes[1].legend()
#     axes[1].set_title("Differences with Smoothed Profile")
#     # Adjust layout for better spacing
#     plt.tight_layout()
