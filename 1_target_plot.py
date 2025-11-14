# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 15:33:25 2024

@author: claire.dussard
"""
import os
import pathlib
from neo.io import AlphaOmegaIO
import matplotlib.pyplot as plt
from functions import get_plot,load_electrode_config,desired_filesV2

os.chdir("../../")
lustre_data_dir = "DATA/"
lustre_path = pathlib.Path(lustre_data_dir)
os.chdir(lustre_path)

noms_patients = os.listdir()[1:]

#lire les donnees
patient_ind = 1
patient_name = noms_patients[patient_ind]

print(patient_name)

#cotes
side = "left"
cote = "/"+side+"/"

lazy = False
types_data = ["LFP","RAW","SPK"]
type_a_afficher = "LFP"
ind_type = types_data.index(type_a_afficher)

targets_df, controls_df = load_electrode_config()

pt_targets_df = targets_df[(targets_df['num_dossier'] == patient_ind+2) & (targets_df['side'] == side)]
pt_controls_df = controls_df[(controls_df['num_dossier'] == patient_ind+2) & (controls_df['side'] == side)]

targets = [(row['electrode'], [row['depth']]) for _, row in pt_targets_df.iterrows()]
controls = [(row['electrode'], [row['depth']]) for _, row in pt_controls_df.iterrows()]



fichiers_patient_target = AlphaOmegaIO(patient_name+cote).read()[0]

fmin = 3
fmax = 80
files_target,desc_target,files_target_filtered = desired_filesV2(fichiers_patient_target,
                                            targets+controls,patient_name,cote,fmin,fmax,ind_type)


filter_data = True
fmin_disp = 6
fmax_disp = 45
if filter_data:
    files = files_target_filtered
else:
    files = files_target
    
print(len(files))
    
# Trier desc_target et récupérer les indices du tri
sorted_indices = sorted(range(len(desc_target)), key=lambda i: desc_target[i][1], reverse=True)
desc_sorted = [desc_target[i] for i in sorted_indices]
files_sorted = [files[i] for i in sorted_indices]

sameYscale = False
get_plot(fmin_disp, fmax_disp, files_sorted, desc_sorted, sameYscale)

#◘===================== SAVE IT ALL ===============

targets_df, controls_df = load_electrode_config()

noms_patients = os.listdir()[1:]

# Lire les données
lazy = False
save = True
types_data = ["LFP", "RAW", "SPK"]
type_a_afficher = "LFP"
ind_type = types_data.index(type_a_afficher)

# Plage de fréquences pour le filtrage et l'affichage
fmin = 2
fmax = 100
fmin_disp = 6
fmax_disp = 50

# Loop sur les patients et les côtés
for patient_ind in range(len(noms_patients)):  # 0 à 4 inclus
    if patient_ind == 5: 
        continue
    else:
        print(patient_ind)
        for side in ["right", "left"]:
            # Définir le nom du patient et le côté
            patient_name = noms_patients[patient_ind]
            print(patient_name)
            cote = "/" + side + "/"
            print(cote)
            # Définir les cibles et contrôles
            pt_targets_df = targets_df[(targets_df['num_dossier'] == patient_ind+2) & (targets_df['side'] == side)]
            pt_controls_df = controls_df[(controls_df['num_dossier'] == patient_ind+2) & (controls_df['side'] == side)]
            print(pt_targets_df)
            print(pt_controls_df)
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
    
            # Choisir les fichiers filtrés ou non
            filter_data = True
            if filter_data:
                files = files_target_filtered
            else:
                files = files_target
    
            # Trier les fichiers et descriptions
            sorted_indices = sorted(range(len(desc_target)), key=lambda i: desc_target[i][1], reverse=True)
            desc_sorted = [desc_target[i] for i in sorted_indices]
            files_sorted = [files[i] for i in sorted_indices]
    
            # Générer et sauvegarder le plot
            filename = f"{patient_name}_{side}.png"
            get_plot(fmin_disp, fmax_disp, files_sorted, desc_sorted, True)
            if save:
                plt.savefig("../img/"+filename)
                plt.close()
            
            filename = f"{patient_name}_{side}_unscaled.png"
            get_plot(fmin_disp, fmax_disp, files_sorted, desc_sorted, False)
            if save:
                plt.savefig("../img/"+filename)
                plt.close() 
    
