# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:43:25 2024

@author: claire.dussard
"""

# General imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Import some NeuroDSP functions to use with MNE
from neurodsp.spectral import compute_spectrum
from specparam import SpectralModel
from specparam.plts.spectra import plot_spectra
from neurodsp.burst import detect_bursts_dual_threshold

# Import NeuroDSP plotting functions
from neurodsp.plts import (plot_time_series, plot_power_spectra,plot_bursts)

from neurodsp.filt import filter_signal

from neurodsp.utils import create_times


fs = 1375.0

    
def load_electrode_config():

    #config_df = pd.read_csv("../elec_targets.csv")#,sep=";")
    config_df = pd.read_excel("../targets.xlsx")
    targets_df = config_df[config_df['type'] == 'target']
    controls_df = config_df[config_df['type'] == 'control']
    return targets_df, controls_df

def desired_filesV2(block,descriptions,patient_name,cote,fmin,fmax,ind_type):
    final_files = []
    files_filtered = []
    descr = []
    for segment in block.segments: #loop sur les enregistrements 
        ini = len(patient_name+cote)+4
        depth = float(segment.file_origin[ini:ini+5])
        analog_signals_type = segment.analogsignals[ind_type]
        for target_electrode, depths in descriptions:
            for k, spiketrain in enumerate(segment.spiketrains):# Loop over each electrode in segment
                electrode_name = spiketrain.name[9:]  # Extract electrode name
                print(electrode_name)
                signal_raw = np.array(analog_signals_type)[:,k]
                if electrode_name == target_electrode and depth in depths:
                    print("trouve")
                    final_files.append(signal_raw )
                    sig_filt = filter_files(signal_raw ,fmin,fmax)
                    files_filtered.append(sig_filt)
                    descr.append((electrode_name,depth))
    return final_files,descr,files_filtered


def filter_files(data_files,fmin,fmax):
    
    filtered_data_files = filter_signal(data_files, fs, 'bandpass', f_range = [fmin,fmax],
                             filter_type='iir',butterworth_order=3,remove_edges=False)
    
    return filtered_data_files
        
def get_plot(fmin,fmax,files,desc,sameYscale):
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
        f_range = (13, 35)

        # Plot time series
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
    
def get_plotV2(fmin, fmax, files, desc, sameYscale,peak_width_limits,min_peak_height,fmax_analyse):
    num_signals = len(files)
    fig, axes = plt.subplots(num_signals, 2, figsize=(14, 5 * num_signals), sharex=False)

    global_min_psd, global_max_psd = float('inf'), float('-inf')
    global_min_sig, global_max_sig = float('inf'), float('-inf')

    for sig in files:
        _, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])
        global_min_psd = min(global_min_psd, np.min(psd_med))
        global_max_psd = max(global_max_psd, np.max(psd_med))
        global_min_sig = min(global_min_sig, np.min(sig))
        global_max_sig = max(global_max_sig, np.max(sig))

    for i, (sig, descr) in enumerate(zip(files, desc)):
        times = create_times(len(sig) / fs, fs)
        if len(times) != len(sig):
            min_length = min(len(times), len(sig))
            times = times[:min_length]
            sig = sig[:min_length]

        freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])
        
        ax = axes[i,0]# if num_signals > 1 else axes  # Handle single subplot case
        ax.plot(freq_med, psd_med)
        print(freq_med[np.argmax(psd_med)])
        ax.plot(freq_med[np.argmax(psd_med)], np.max(psd_med), '.r', ms=12)
        
        fm = SpectralModel(aperiodic_mode="knee", min_peak_height=min_peak_height, peak_width_limits=peak_width_limits, max_n_peaks=6)
        fm.fit(freq_med, psd_med, (fmin, fmax_analyse))
        fm.plot(ax=axes[i, 1])
        #plot_spectra(fm.freqs, fm.modeled_spectrum_, ax=axes[i, 1], label='Full Model', color='red')
        axes[i, 1].set_title(f"{descr} | Error: {fm.get_params('error'):.2f}")

        if sameYscale:
            axes[i, 0].set_ylim(global_min_sig, global_max_sig)
            axes[i, 1].set_ylim(global_min_psd, global_max_psd)

    plt.tight_layout()
    plt.show()

    
# def desired_files(list_files,descriptions,patient_name,cote,fmin,fmax,ind_type):
#     final_files = []
#     files_filtered = []
#     descr = []
#     n_segments = len(list_files._segments)
#     for i in range(n_segments): #loop sur les enregistrements 
#         seg_i = list_files.segments[i]
#         ini = len(patient_name+cote)+4
#         profondeur = float(seg_i.file_origin[ini:ini+5])
#         data_record_i = seg_i.analogsignals
#         data_record_i_type_j = data_record_i[ind_type]
#         for target in descriptions:
#             elec = target[0]
#             ls_depths = target[1]
#             if (profondeur in ls_depths):
#                 for k in range(len(seg_i.spiketrains)):#loop sur les electrodes
#                     data_record_i_type_j_elec_k = np.array(data_record_i_type_j)[:,k]
#                     print(seg_i.spiketrains[k].name[9:])
#                     if (seg_i.spiketrains[k].name[9:]==elec):
#                         print("found "+elec+str(profondeur))
#                         final_files.append(data_record_i_type_j_elec_k)
 
#                         sig_filt = filter_files(data_record_i_type_j_elec_k,fmin,fmax)
#                         files_filtered.append(sig_filt)
#                         descr.append((elec,profondeur))
#     return final_files,descr,files_filtered

