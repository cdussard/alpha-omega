# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:43:25 2024

@author: claire.dussard
"""

# General imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors

# Import MNE, as well as the MNE sample dataset
import mne
from mne import io
from mne.datasets import sample
from mne.viz import plot_topomap

# Import some NeuroDSP functions to use with MNE
from neurodsp.spectral import compute_spectrum, trim_spectrum
from neurodsp.burst import detect_bursts_dual_threshold
from neurodsp.rhythm import compute_lagged_coherence

# Import NeuroDSP plotting functions
from neurodsp.plts import (plot_time_series, plot_power_spectra,
                           plot_bursts, plot_lagged_coherence)

from neurodsp.filt import filter_signal
from neurodsp.spectral import compute_spectrum, rotate_powerlaw

# Import utilities for loading and plotting data
from neurodsp.utils import create_times
from neurodsp.utils.download import load_ndsp_data
from neurodsp.plts.spectral import plot_power_spectra
from neurodsp.plts.time_series import plot_time_series
# Import some NeuroDSP functions to use with MNE
from neurodsp.spectral import compute_spectrum, trim_spectrum
from neurodsp.burst import detect_bursts_dual_threshold
from neurodsp.rhythm import compute_lagged_coherence

# Import NeuroDSP plotting functions
from neurodsp.plts import (plot_time_series, plot_power_spectra,
                           plot_bursts, plot_lagged_coherence)


fs = 1375.0


def desired_files(list_files,descriptions,ctrl,patient_name,cote,fmin,fmax):
    ind_type = 0
    final_files = []
    files_filtered = []
    descr = []
    n_segments = len(list_files._segments)
    for i in range(n_segments): #loop sur les enregistrements 
        seg_i = list_files.segments[i]
        ini = len(patient_name+cote)+4
        if ctrl:
            ini = ini + 7
        profondeur = float(seg_i.file_origin[ini:ini+5])
        data_record_i = seg_i.analogsignals
        data_record_i_type_j = data_record_i[ind_type]
        for target in descriptions:
            elec = target[0]
            ls_depths = target[1]
            if (profondeur in ls_depths):
                for k in range(len(seg_i.spiketrains)):#loop sur les electrodes
                    data_record_i_type_j_elec_k = np.array(data_record_i_type_j)[:,k]
                    print(seg_i.spiketrains[k].name[9:])
                    if (seg_i.spiketrains[k].name[9:]==elec):
                        print("FOUND "+elec+str(profondeur))
                        
                        final_files.append(data_record_i_type_j_elec_k)
 

                        sig_filt2 = filter_signal(data_record_i_type_j_elec_k, 1375, 'bandpass', f_range = [fmin,fmax],#2,80
                                                 filter_type='iir',butterworth_order=3,remove_edges=False)
                        files_filtered.append(sig_filt2)
                        descr.append((elec,profondeur))
    return final_files,descr,files_filtered
        
def get_plot(fmin,fmax,files,desc,sameYscale):
   
    # Number of signals
    num_signals = len(files)
    
    # Create a figure with subplots arranged in two columns
    fig, axes = plt.subplots(num_signals, 2, figsize=(14, 5 * num_signals), sharex=False)
    
    global_min = float('inf')
    global_max = float('-inf')
    
    for sig in files:
        # Compute the power spectrum
        _, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])
        global_min = min(global_min, np.min(psd_med))
        global_max = max(global_max, np.max(psd_med))
        global_min_sig = min(global_min, np.min(sig))
        global_max_sig = max(global_max, np.max(sig))
        
    # Loop over signals and their descriptions
    for i, (sig, descr) in enumerate(zip(files, desc)):
        print(sig)
        print(descr)
        times = create_times(len(sig) / fs, fs)
        # Check if the lengths of `times` and `sig` match
        if len(times) != len(sig):
            min_length = min(len(times), len(sig))
            times = times[:min_length]
            sig = sig[:min_length]
        # Plot the loaded signal
        amp_dual_thresh = (1, 2)
        f_range = (8, 12)

        # Detect bursts using dual threshold algorithm
        bursting = detect_bursts_dual_threshold(sig, fs, amp_dual_thresh, f_range)
        plot_bursts(times, sig, bursting, labels=[descr],ax=axes[i,0])
        # Compute the power spectrum
        freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])
    
        # Plot on the appropriate subplot
        ax = axes[i,1]# if num_signals > 1 else axes  # Handle single subplot case
        ax.plot(freq_med, psd_med)
        print(freq_med[np.argmax(psd_med)])
        ax.plot(freq_med[np.argmax(psd_med)], np.max(psd_med), '.r', ms=12)
        ax.set_title(descr)  # Title each subplot
        if sameYscale:
            axes[i,0].set_ylim(global_min_sig, global_max_sig)
            ax.set_ylim(global_min, global_max)
    
    # Adjust layout
    plt.tight_layout()
    plt.show()