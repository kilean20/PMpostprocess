from .utils import *
from .profilemonitor_plot_util import profilemonitor_plot as pmplt
import os
import re


def list_raw_files(data_path):
    path = os.path.join(data_path, 'data')
    pattern = re.compile(r".*_PM_D\d{4}_\d{8}_\d{6}\.dat$")
    files = [f for f in os.listdir(path) if pattern.match(f)]
    return files


def read_raw_file(fname):
    # Retrieve parameters and raw data
    params = pmplt.get_param(fname)
    raw = pmplt.get_rawdata(fname, params)
    
    data = []
    
    for i in range(1, 4):  # Process for indices 1, 2, 3
        # Extract and sort positions and signals
        pos = np.array(raw[f'a1p{i}'].split()[4:], dtype=float)
        sig = np.array(raw[f'a1s{i}'].split()[4:], dtype=float)
        idx = np.argsort(pos)
        pos = pos[idx]
        sig = sig[idx]
        
        # Identify largest gaps in positions
        gaps = pos[1:] - pos[:-1]
        largest_gaps = np.sort(np.argsort(gaps)[-2:])
        
        # Find the domain containing the peak signal
        peak = np.argmax(sig)
        
        if largest_gaps[1] < peak:  # Peak is after the second largest gap
            pos = pos[largest_gaps[1] + 1:]
            sig = sig[largest_gaps[1] + 1:]
        elif largest_gaps[0] < peak:  # Peak is between the two largest gaps
            pos = pos[largest_gaps[0] + 1:largest_gaps[1]]
            sig = sig[largest_gaps[0] + 1:largest_gaps[1]]
        else:  # Peak is before the first largest gap
            pos = pos[:largest_gaps[0]]
            sig = sig[:largest_gaps[0]]
        # Store the filtered data
        data.append(np.stack((pos, sig)))
    return data


def read_all_raw_data(data_path):
    files = list_raw_files(data_path)
    all_data = {}
    for f in files:
        raw_data = read_raw_file(os.path.join(data_path, 'data', f))
        all_data[f] = raw_data   
    return all_data
