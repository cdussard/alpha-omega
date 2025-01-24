# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 19:42:06 2024

@author: claire.dussard
"""

import numpy as np

files = files_ctrl
desc = descr_ctrl

files = files_target_filtered
desc = desc_target

fmin = 4
fmax = 100


# Number of signals
num_signals = len(files)

global_min = float('inf')
global_max = float('-inf')

for sig in files:
    # Compute the power spectrum
    _, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 3, f_range=[fmin, 40])
    global_min = min(global_min, np.min(psd_med))
    global_max = max(global_max, np.max(psd_med))
    
    
fig, axes = plt.subplots(1,2,figsize=(14, 5 ), sharex=False)


el_c = np.where([desc_[0]=='Central' for desc_ in desc])[0]
files_c = [files[i] for i in el_c]
el_p = np.where([desc_[0]=='Posterior' for desc_ in desc])[0]
files_p = [files[i] for i in el_p]


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

    freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 1.2, f_range=[fmin, fmax])
    ind_beta = np.where(freq_med==35)[0][0]
    # Plot on the appropriate subplot
    ax = axes# if num_signals > 1 else axes  # Handle single subplot case
    #ax.plot(freq_med, psd_med)
    ax[0].plot(freq_med_, psd_med[0:ind_beta], label=descr)
    
    #ax.set_ylim(global_min, global_max)

axes[0].legend(title="Uncorrected")
axes[1].legend(title="Corrected")

# Adjust layout
plt.tight_layout()
plt.show()