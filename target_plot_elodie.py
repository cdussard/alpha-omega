# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 15:33:25 2024

@author: claire.dussard
"""
import os
import pathlib
from neo.io import AlphaOmegaIO
import matplotlib.pyplot as plt
from functions import desired_files,get_plot

os.chdir("../../")
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

fmin = 2
fmax = 40
files_target,desc_target,files_target_filtered = desired_files(fichiers_patient_target,
                                            [targets[patient_ind]]+[controls[patient_ind]],False,patient_name,cote,fmin,fmax)


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


get_plot(fmin_disp, fmax_disp, files_sorted, desc_sorted, True)


#◘===================== SAVE IT ALL ===============
noms_patients = os.listdir()[1:]

# Lire les données
lazy = False
types_data = ["LFP", "RAW", "SPK"]
type_a_afficher = "LFP"
ind_type = types_data.index(type_a_afficher)

# Plage de fréquences pour le filtrage et l'affichage
fmin = 4
fmax = 40
fmin_disp = 6
fmax_disp = 50

# Loop sur les patients et les côtés
for patient_ind in range(5):  # 0 à 4 inclus
    for side in ["right", "left"]:
        # Définir le nom du patient et le côté
        patient_name = noms_patients[patient_ind]
        cote = "/" + side + "/"

        # Définir les cibles et contrôles
        if side == "right":
            targets = [("Posterior", [0.0]), ("Central", [1.0]), ("Medial", [0.0]), ("Central", [1.0]), ("Medial", [2.0])]
            controls = [("Posterior", [4, -1]), ("Central", [3.4, 0]), ("Medial", [3.4, -1.0]), ("Central", [3.0, -0.6]), ("Medial", [4.4, 0.0])]
        elif side == "left":
            targets = [("Central", [0.0]), ("Central", [1.0]), ("Central", [0.0]), ("Central", [0.0]), ("Lateral", [0.0])]
            controls = [("Central", [4, -1]), ("Central", [4.0, -0.4]), ("Central", [3.0, -0.6]), ("Central", [3.4, -0.6]), ("Lateral", [4.0, -0.2])]

        # Lire les fichiers pour le patient et le côté
        fichiers_patient_target = AlphaOmegaIO(patient_name + cote).read()[0]
        files_target, desc_target, files_target_filtered = desired_files(
            fichiers_patient_target,
            [targets[patient_ind]] + [controls[patient_ind]],
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
        filename = f"{patient_name}_{side}_hfo.png"
        get_plot(fmin_disp, fmax_disp, files_sorted, desc_sorted, True)
        plt.savefig(filename)
        plt.close()  # Fermer la figure pour libérer la mémoire
        
        filename = f"{patient_name}_{side}_unscaled_hfo.png"
        get_plot(fmin_disp, fmax_disp, files_sorted, desc_sorted, False)
        plt.savefig(filename)
        plt.close()  # Fermer la figure pour libérer la mémoire

