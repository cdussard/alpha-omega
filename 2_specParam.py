# -*- coding: utf-8 -*-
"""
Created on Thu Jan 9 15:39:57 2025

@author: claire.dussard
"""


import os
import pathlib
from neo.io import AlphaOmegaIO
from functions import *
import pandas as pd
from neurodsp.spectral import compute_spectrum
from specparam import SpectralModel
from specparam.plts.spectra import plot_spectra
from functions import load_electrode_config,desired_filesV2,get_plot,get_plotV2


os.chdir("../../")
lustre_data_dir = "DATA"
lustre_path = pathlib.Path(lustre_data_dir)
os.chdir(lustre_path)

noms_patients = os.listdir()[1:]


fs = 1375
fmax = 400
fmin = 3
fmax_analyse = 50
fmax_disp = 40

#lire les donnees
patient_ind = 1
patient_name = noms_patients[patient_ind]
#cotes
side = "left"
cote = "/"+side+"/"

lazy = False
types_data = ["LFP","RAW","SPK"]
type_a_afficher = "LFP"
ind_type = types_data.index(type_a_afficher)

targets_df, controls_df = load_electrode_config()

# pt_targets_df = targets_df[(targets_df['patient'] == patient_ind) & (targets_df['side'] == side)]
# pt_controls_df = controls_df[(controls_df['patient'] == patient_ind) & (controls_df['side'] == side)]
pt_targets_df = targets_df[(targets_df['num_dossier'] == patient_ind+2) & (targets_df['side'] == side)]
pt_controls_df = controls_df[(controls_df['num_dossier'] == patient_ind+2) & (controls_df['side'] == side)]
print(pt_targets_df)

targets = [(row['electrode'], [row['depth']]) for _, row in pt_targets_df.iterrows()]
controls = [(row['electrode'], [row['depth']]) for _, row in pt_controls_df.iterrows()]
print(targets)
print(controls)

fichiers_patient_target = AlphaOmegaIO(patient_name+cote).read()[0]


files_target,desc_target,files_target_filtered = desired_filesV2(fichiers_patient_target,
                                            targets+controls,patient_name,cote,fmin,fmax,ind_type)


#example de specparam
for sig in files_target:
    freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='mean', nperseg=fs * 2, f_range=[fmin, fmax_analyse])
    fm = SpectralModel(aperiodic_mode="knee",peak_width_limits=[2.5, 20],max_n_peaks=6,min_peak_height = 0.6)
    fm.fit(freq_med, psd_med, (fmin,fmax_analyse))
    fm.print_results()
    fm.plot()
    

#spectre 
sameYscale = False
filter_data = True
fmax_disp = 45
sorted_indices = sorted(range(len(desc_target)), key=lambda i: desc_target[i][1], reverse=True)
desc_sorted = [desc_target[i] for i in sorted_indices]
files_sorted = [files_target[i] for i in sorted_indices]


peak_width_limits = [1,20]
min_peak_height = 0.6
get_plotV2(fmin, fmax_disp, files_target, desc_target, sameYscale,peak_width_limits,min_peak_height,fmax_analyse)


# =========== tous les patients==============
fmin = 3
for patient_ind in range(5):  # 0 à 4 inclus
    for side in ["right", "left"]:
        # Définir le nom du patient et le côté
        patient_name = noms_patients[patient_ind]
        cote = "/" + side + "/"

        # Définir les cibles et contrôles
        pt_targets_df = targets_df[(targets_df['num_dossier'] == patient_ind+2) & (targets_df['side'] == side)]
        pt_controls_df = controls_df[(controls_df['num_dossier'] == patient_ind+2) & (controls_df['side'] == side)]

        targets = [(row['electrode'], [row['depth']]) for _, row in pt_targets_df.iterrows()]
        controls = [(row['electrode'], [row['depth']]) for _, row in pt_controls_df.iterrows()]

        # Lire les fichiers pour le patient et le côté
        fichiers_patient_target = AlphaOmegaIO(patient_name + cote).read()[0]
        
        files_target, desc_target, files_target_filtered = desired_filesV2(
            fichiers_patient_target,
            targets + controls,
            patient_name,
            cote,
            fmin,
            fmax,
            ind_type
        )

        sorted_indices = sorted(range(len(desc_target)), key=lambda i: desc_target[i][1], reverse=True)
        desc_sorted = [desc_target[i] for i in sorted_indices]
        files_sorted = [files_target[i] for i in sorted_indices]
        get_plotV2(fmin, fmax_disp, files_sorted, desc_sorted, sameYscale,peak_width_limits,min_peak_height,fmax_analyse)


#==============EXTRACT THE DATA ===============
# Frequency range for filtering
# Loop over patients and sides
results = []
fs = 1375  # Sampling frequency
for patient_ind in range(len(noms_patients)):  # 0 à 4 inclus
    if patient_ind ==5:
        continue
    for side in ["right", "left"]:
        # Définir le nom du patient et le côté
        patient_name = noms_patients[patient_ind]
        cote = "/" + side + "/"

        # Définir les cibles et contrôles
        pt_targets_df = targets_df[(targets_df['num_dossier'] == patient_ind+2) & (targets_df['side'] == side)]
        pt_controls_df = controls_df[(controls_df['num_dossier'] == patient_ind+2) & (controls_df['side'] == side)]

        targets = [(row['electrode'], [row['depth']]) for _, row in pt_targets_df.iterrows()]
        controls = [(row['electrode'], [row['depth']]) for _, row in pt_controls_df.iterrows()]

        # Lire les fichiers pour le patient et le côté
        fichiers_patient_target = AlphaOmegaIO(patient_name + cote).read()[0]
        
        files_target, desc_target, files_target_filtered = desired_filesV2(
            fichiers_patient_target,
            targets + controls,
            patient_name,
            cote,
            fmin,
            fmax,
            ind_type
        )

        files = files_target

        # Loop through each signal
        for sig, desc in zip(files, desc_target):
            electrode, depth = desc
            # Determine if the electrode is a target or control
            if any(elec == electrode and depth in depths for elec, depths in targets):
                elec_type = "target"
            elif any(elec == electrode and depth in depths for elec, depths in controls):
                elec_type = "control"
            else:
                elec_type = "unknown"

            # Determine the region
            if elec_type == "target":
                region = "STN"
            elif elec_type == "control" and depth > 0:
                region = "above"
            elif elec_type == "control" and depth <= 0:
                region = "SNR"
            else:
                region = "unknown"

            # Compute power
            freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=fs * 2, f_range=[fmin, fmax])
            
            fm = SpectralModel(aperiodic_mode="knee",peak_width_limits=peak_width_limits,min_peak_height=min_peak_height, max_n_peaks=6)
            fm.fit(freq_med, psd_med, (fmin,fmax_analyse))
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
                    "type": elec_type,  # type (target or control)
                    "region": region,  # region (STN, above, SNR)
                })


# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Display the DataFrame
print(results_df)

# Optionally save to a CSV file
results_df.to_csv("../csv_files/specParam_results.csv", index=False)


#examples de visu
import matplotlib.pyplot as plt
import seaborn as sns
df = results_df
df_sub40 = df[df['peak_freq'] < 40]

# Boxplot of Peak Frequency < 40Hz by Patient and Region
plt.figure(figsize=(14, 6))
sns.boxplot(data=df_sub40, x='patient_ind', y='peak_freq', hue='region', palette='Set2')
sns.stripplot(data=df_sub40, x='patient_ind', y='peak_freq', hue='region', dodge=True, color='black', alpha=0.4, jitter=True)

plt.title("Peak Frequency < 40 Hz by Patient and Region")
plt.ylabel("Peak Frequency (Hz)")
plt.xlabel("Patient Index")
plt.legend(title='Region', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()



plt.figure(figsize=(10, 6))
sns.violinplot(data=df, x='region', y='exponent', inner=None, palette='Pastel1', cut=0, alpha=0.6)
sns.stripplot(data=df, x='region', y='exponent', color='black', alpha=0.5, jitter=True)

# Median lines
medians_exp = df.groupby('region')['exponent'].median()
for i, region in enumerate(medians_exp.index):
    plt.plot([i - 0.2, i + 0.2], [medians_exp[region]]*2, color='red', linewidth=2)

plt.title("Aperiodic Exponent Distribution by Brain Region\n(Violin + Transparent Points + Median Line)")
plt.ylabel("Aperiodic Exponent")
plt.xlabel("Region")
plt.tight_layout()
plt.show()