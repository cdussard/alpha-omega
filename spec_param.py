# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:48:52 2025

@author: claire.dussard
"""

import os
import pathlib
from neo.io import AlphaOmegaIO
from neo.core import Block, Segment, Group, AnalogSignal
import matplotlib.pyplot as plt
import numpy as np
from functions import *
import pandas as pd
import specparam
from neurodsp.spectral import compute_spectrum, rotate_powerlaw
from specparam import SpectralModel
# Import the FOOOF object
from fooof import FOOOF

os.chdir("../../")
lustre_data_dir = "DATA/"
lustre_path = pathlib.Path(lustre_data_dir)
os.chdir(lustre_path)

noms_patients = os.listdir()[1:]

#lire les donnees
patient_ind = 3
patient_name = noms_patients[patient_ind]

print(patient_name)

#cotes
side = "right"
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


files_target,desc_target,files_target_filtered = desired_files(fichiers_patient_target,targets,False,patient_name,cote,fmin,500)


# Loop over signals and their descriptions
fs = 1375

fmin = 1
fmax = 400
fmax_analyse = 350

for sig in files_target:

    # Compute the power spectrum
    freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='mean', nperseg=fs * 2, f_range=[fmin, fmax])
    

    fm = FOOOF(aperiodic_mode="knee",peak_width_limits=[2.5, 50],max_n_peaks=6)
    fm.fit(freq_med, psd_med, (fmin,fmax_analyse))
    print(freq_med[1]-freq_med[0])#frequency resolution
    fm.print_results()
    fm.plot()
    # Obtenir les résultats
    ap_params, peak_params, r_squared, fit_error, gauss_params = fm.get_results()

    
    
#% =================== FOR ALL PATIENTS ==============
# Frequency range for filtering
fmin, fmax = 1, 80  # Example values (update as needed)

# Loop over patients and sides
results = []
fs = 1375  # Sampling frequency

for patient_ind in range(len(noms_patients)):  # Patients 0 to 4
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
            cote=cote,fmin=fmin,fmax=150
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
            elif elec_type == "control" and depth <= 0:
                region = "SNR"
            else:
                region = "unknown"

            fmax = 70
            # Compute the power spectrum
            freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='mean', nperseg=fs * 2, f_range=[fmin, fmax])
            
        
            fm = FOOOF(aperiodic_mode="knee",peak_width_limits=[2.5, 20],max_n_peaks=6)
            fmax_analyse = 60
            fm.fit(freq_med, psd_med, (fmin,fmax_analyse))
            # Obtenir les résultats
            ap_params, peak_params, r_squared, fit_error, gauss_params = fm.get_results()

            # Append results for each peak to the list
            for peak_param in peak_params:
                results.append({
                    "side": side,
                    "patient_ind": patient_ind,
                    "peak_freq": peak_param[0],
                    "peak_height": peak_param[1],
                    'peak_width': peak_param[2],
                    "electrode": electrode,
                    "depth": depth,
                    'aperiodic_offset': ap_params[0],  # Offset aperiodique
                    'knee': ap_params[1],              # Paramètre knee
                    'exponent': ap_params[2],          # Exposant
                    'r_squared': r_squared,            # R-squared
                    'fit_error': fit_error ,
                    'gauss_params': gauss_params,      # Paramètres Gaussien (s'il y en a)
                    "type": elec_type,  # Add the type (target or control)
                    "region": region,  # Add the region (STN, above, SNR)
                })


# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Display the DataFrame
print(results_df)

# Optionally save to a CSV file
results_df.to_csv("../specParam_1_80Hz.csv", index=False)



    