# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 19:57:29 2025

@author: claire.dussard
"""


import numpy as np
from scipy.signal import find_peaks, peak_widths

import os
import pathlib
from neo.io import AlphaOmegaIO
from neo.core import Block, Segment, Group, AnalogSignal
import matplotlib.pyplot as plt
import numpy as np
from functions import *

os.chdir("./")
lustre_data_dir = "DATA/"
lustre_path = pathlib.Path(lustre_data_dir)
os.chdir(lustre_path)

noms_patients = os.listdir()[1:]

#lire les donnees
patient_ind = 4
patient_name = noms_patients[patient_ind]

print(patient_name)

#cotes
side = "left"
cote = "/"+side+"/"

lazy = False
types_data = ["LFP","RAW","SPK"]
type_a_afficher = "LFP"
ind_type = types_data.index(type_a_afficher)

#A TERME LES LIRE DANS TABLEUR
if side =="right":
    targets = [("Posterior",[0.0]),("Central",[1.0]),("Medial",[0.0]),("Central",[1.0]),("Medial",[2.0])]
    controls = [("Posterior",[4,-1]),("Central",[3.4,0]),("Medial",[3.4,-1.0]),("Central",[3.0,-0.6]),("Medial",[4.4,0.0])]
elif side =="left":
    targets = [("Central",[0.0]),("Central",[1.0]),("Central",[0.0]),("Central",[0.0]),("Lateral",[0.0])]
    controls = [("Central",[4,-1]),("Central",[4.0,-0.4]),("Central",[3.0,-0.6]),("Central",[3.4,-0.6]),("Lateral",[4.0,-0.2])]#premier = debut STN, deuxieme = SNr, obligé arrondir au nb pair 
    
print(targets[patient_ind])
print(controls[patient_ind])

fichiers_patient_target = AlphaOmegaIO(patient_name+cote).read()[0]


files_target,desc_target,files_target_filtered = desired_files(fichiers_patient_target,targets,False,patient_name,cote)


# Loop over signals and their descriptions
fs = 1375

sig = files_target[0]

num_signals = 2



# Compute the power spectrum
_, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])

    
descr = sorted(descr_ctrl, key=lambda x: x[1], reverse=True)[0]


freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])

peaks, properties = find_peaks(psd_med, height=0.01,distance=7,width=1)
print(freq_med[peaks])


noms_patients = os.listdir()[1:]

import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths

# Frequency range for filtering
fmin, fmax = 2, 80  # Example values (update as needed)

# Loop over patients and sides
results = []
fs = 1375  # Sampling frequency
for patient_ind in range(5):  # Patients 0 to 4
    for side in ["right", "left"]:
        # Define patient name and side
        patient_name = noms_patients[patient_ind]
        cote = "/" + side + "/"

        # Define targets and controls
        if side == "right":
            targets = [("Posterior", [0.0]), ("Central", [1.0]), ("Medial", [0.0]), ("Central", [1.0]), ("Medial", [2.0])]
            controls = [("Posterior", [4, -1]), ("Central", [3.4, 0]), ("Medial", [3.4, -1.0]), ("Central", [3.0, -0.6]), ("Medial", [4.4, 0.0])]
        elif side == "left":
            targets = [("Central", [0.0]), ("Central", [1.0]), ("Central", [0.0]), ("Central", [0.0]), ("Lateral", [0.0])]
            controls = [("Central", [4, -1]), ("Central", [4.0, -0.4]), ("Central", [3.0, -0.6]), ("Central", [3.4, -0.6]), ("Lateral", [4.0, -0.2])]

        # Read files for the patient and side
        fichiers_patient_target = AlphaOmegaIO(patient_name + cote).read()[0]

        # Use the desired_files function to extract signals
        files_target, desc_target, files_target_filtered = desired_files(
            fichiers_patient_target,
            [targets[patient_ind]] + [controls[patient_ind]],
            ctrl=False,
            patient_name=patient_name,
            cote=cote
        )

        # Choose whether to use filtered or raw signals
        filter_data = True
        if filter_data:
            files = files_target_filtered
        else:
            files = files_target

        # Loop through each signal
        for sig, desc in zip(files, desc_target):
            # Extract electrode and depth from description
            electrode = desc[0]
            depth = desc[1]

            # Determine if the electrode is a target or control
            if any(elec == electrode and depth in depths for elec, depths in targets):
                elec_type = "target"
            elif any(elec == electrode and depth in depths for elec, depths in controls):
                elec_type = "control"
            else:
                elec_type = "unknown"  # In case it's neither target nor control

            # Determine the region based on type and depth
            if elec_type == "target":
                region = "STN"
            elif elec_type == "control" and depth > 0:
                region = "above"
            elif elec_type == "control" and depth < 0:
                region = "SNR"
            else:
                region = "unknown"

            # Compute the power spectrum
            freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])

            # Find peaks in the PSD
            peaks, properties = find_peaks(psd_med, height=0.01, distance=7, width=1)

            # Calculate peak widths
            widths, width_heights, left_ips, right_ips = peak_widths(psd_med, peaks, rel_height=0.5)

            # Extract peak frequencies, prominence, and other properties
            peak_frequencies = freq_med[peaks]
            prominences = properties['prominences']

            # Append results for each peak to the list
            for i in range(len(peaks)):
                results.append({
                    "side": side,
                    "patient_ind": patient_ind,
                    "peak_freq": peak_frequencies[i],
                    "prominence": prominences[i],
                    "width": widths[i],
                    "width_heights": width_heights[i],
                    "electrode": electrode,
                    "depth": depth,
                    "type": elec_type,  # Add the type (target or control)
                    "region": region,  # Add the region (STN, above, SNR)
                })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Display the DataFrame
print(results_df)

# Optionally save to a CSV file
results_df.to_csv("peaks_results.csv", index=False)
