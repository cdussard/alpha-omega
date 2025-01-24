# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:08:14 2024

@author: claire.dussard
"""

import os
import pathlib
from neo.io import AlphaOmegaIO
from neo.core import Block, Segment, Group, AnalogSignal
import matplotlib.pyplot as plt
import numpy as np



os.chdir("./")
lustre_data_dir = "DATA/"
lustre_path = pathlib.Path(lustre_data_dir)
os.chdir(lustre_path)

noms_patients = os.listdir()

#lire les donnees
patient_ind = 1
patient_name = noms_patients[patient_ind]

print(patient_name)
#profondeurs a trouver

#cotes
side = "right"
cote = "/"+side+"/"

lazy = False
types_data = ["LFP","RAW","SPK"]

fichiers_patient = AlphaOmegaIO(patient_name+cote).read()[0]#lsx_files=["iT1D10.000F0001.lsx"]

#tu as n_segments = nombre de fichiers + 1 
n_segments = len(fichiers_patient._segments)

for i in range(n_segments): #loop sur les enregistrements 
    seg_i = fichiers_patient.segments[i]
    print(seg_i.name)
    ini = len(patient_name+cote)+4
    profondeur = seg_i.file_origin[ini:ini+5]
    print(profondeur)
    data_record_i = seg_i.analogsignals
    for j in range(3):# 3 types : RAW LFP ET SPIKE
        data_record_i_type_j = data_record_i[j]
        print(data_record_i_type_j.name)
        for k in range(len(seg_i.spiketrains)):#loop sur les electrodes
            data_record_i_type_j_elec_k = np.array(data_record_i_type_j)[:,k]
            print(seg_i.spiketrains[k].name)
            plt.plot(data_record_i_type_j.times, data_record_i_type_j)
            
            
            
            
type_a_afficher = "LFP"
ind_type = types_data.index(type_a_afficher)
#avec les subplots

# Création de la figure et des subplots
fig, axs = plt.figure(figsize=(15, 10))  # 3 types, n_segments colonnes

# Boucle sur les segments
for i in range(n_segments):
    seg_i = fichiers_patient.segments[i]
    print(seg_i.name)
    
    ini = len(patient_name + cote) + 4
    profondeur = seg_i.file_origin[ini:ini + 5]
    print(profondeur)
    
    data_record_i = seg_i.analogsignals
    #Prend le bon type (RAW, LFP, SPIKE)
    data_record_i_type = data_record_i[ind_type]
    print(data_record_i_type.name)
    
    for j in range(len(seg_i.spiketrains)):        # Boucle sur les électrodes
        ax = axs[j, i]
        data_record_i_type_elec_k = np.array(data_record_i_type_j)[:,j]#MARCHE PAS 
        print(seg_i.spiketrains[j].name)
        
        # Plot sur le subplot correspondant
        ax.plot(data_record_i_type.times, data_record_i_type_elec_k)
        ax.set_title(f"{profondeur} mm, {type_a_afficher},{seg_i.spiketrains[j].name[8:]}")

  

#tu as n_segments = nombre de fichiers + 1 
n_segments = len(fichiers_patient._segments)

for i in range(n_segments): #loop sur les enregistrements 
    seg_i = fichiers_patient.segments[i]
    print(seg_i.name)
    ini = len(patient_name+cote)+4
    profondeur = seg_i.file_origin[ini:ini+5]
    print(profondeur)
    data_record_i = seg_i.analogsignals
    for j in range(3):# 3 types : RAW LFP ET SPIKE
        data_record_i_type_j = data_record_i[j]
        print(data_record_i_type_j.name)
        for k in range(len(seg_i.spiketrains)):#loop sur les electrodes
            data_record_i_type_j_elec_k = np.array(data_record_i_type_j)[:,k]
            print(seg_i.spiketrains[k].name)
            plt.plot(data_record_i_type_j.times, data_record_i_type_j)
            




data_record_i_type_j = data_record_i[0]#LFP
w = data_record_i_type_j
data = np.array(w)
fs = int(data_record_i_type_j.sampling_rate)

for i in range(3):#un spectre par electrode
 # Convertit l'objet AnalogSignal en tableau Numpy
    column_2 = data[:, i]
    sig = column_2
    
    # Set sampling rate, and create a times vector for plotting
    times = create_times(len(column_2)/fs, fs)
    # Plot the loaded signal
    plot_time_series(times, column_2)
    
    # Median of spectrogram ("median Welch")
    freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs*2)

    plot_power_spectra([ freq_med[:200]],
                       [ psd_med[:200]],
                       ['Median Welch'])
    

    freqs, powers = trim_spectrum(freq_med, psd_med, [3, 30])
    # Check where the peak power is
    peak_cf = freqs[np.argmax(powers)]
    print(peak_cf)
    plot_power_spectra(freqs, powers)
    plt.plot(freqs[np.argmax(powers)], np.max(powers), '.r', ms=12)
    
    # Burst settings
    amp_dual_thresh = (1., 1.5)
    f_range = (peak_cf-2, peak_cf+2)
    # Detect bursts of high amplitude oscillations in the extracted signal
    bursting = detect_bursts_dual_threshold(sig, fs, amp_dual_thresh, f_range)
    # Plot original signal and burst activity
    plot_bursts(times, sig, bursting, labels=['Raw Data', 'Detected Bursts'])



