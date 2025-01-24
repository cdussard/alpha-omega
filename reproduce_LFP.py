# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 12:31:03 2024

@author: claire.dussard
"""

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

noms_patients = os.listdir()

#lire les donnees
patient_ind = 1
patient_name = noms_patients[patient_ind]

print(patient_name)
#droit

#cotes
side = "right"
cote = "/"+side+"/"

lazy = False
types_data = ["LFP","RAW","SPK"]
type_a_afficher = "LFP"
ind_type = types_data.index(type_a_afficher)

targets = [("Central",[0.8,0.4,0]),("Posterior",[-0.2,-0.4,-0.8])]

targets = [("Central",[0.8,0.4,0]),("Posterior",[0.8,0.4,0])]
controls = [("Lateral",[0.4,5.0,2.0])]

#VD pat
targets = [("Central",[1.0,0.8,1.2])]
controls = [("Central",[4.4,5.0,0.0,-0.2])]


#-controls_outside = [("Central",[10.0,5.0,8.7]),("Posterior",[10.0,5.0,8.7])]

fichiers_patient_target = AlphaOmegaIO(patient_name+cote).read()[0]

#fichiers_patient_controls = AlphaOmegaIO(patient_name+cote+"/check2/").read()[0]#on ne sait pas pk mais les RT2 ne passent pas

#tu as n_segments = nombre de fichiers + 1 
n_segments = len(fichiers_patient_target._segments)


#======================get all files==============================

files_target,desc_target,files_target_filtered = desired_files(fichiers_patient_target,targets,False,patient_name,cote)

#
files_ctrl,descr_ctrl,files_ctrl_filtered = desired_files(fichiers_patient_target,controls,False,patient_name,cote)

#files_ctrl,descr_ctrl,files_ctrl_filtered = desired_files(fichiers_patient_controls,controls_outside,True,patient_name,cote)

#=====================GET THE PLOTS ================

fmin = 6.5
fmax = 45
get_plot(fmin,fmax,files_target,desc_target,False)

# desc_sorted = sorted(desc_target, key=lambda x: x[0],reverse=False)
# desc_ctrl_sorted = sorted(descr_ctrl, key=lambda x: x[0],reverse=False)

desc_sorted = sorted(desc_target, key=lambda x: x[1], reverse=True)
desc_ctrl_sorted = sorted(descr_ctrl, key=lambda x: x[1], reverse=True)

#get_plot(fmin,fmax,files_target_filtered,desc_sorted,False)
get_plot(fmin,fmax,files_target_filtered,desc_sorted,True)


#get_plot(fmin,fmax,files_ctrl_filtered,desc_ctrl_sorted,False)
get_plot(fmin,fmax,files_ctrl_filtered,desc_ctrl_sorted,True)

#==========================================GET A TFR??==========================================
import mne

files = files_target_filtered

desc_ = desc_target
files = files_target

desc_ = descr_ctrl
files = files_ctrl

ind = 0
str_ = desc_[ind][0]+str(desc_[0][1])

# Assuming `files_target[ind]` contains your data and you know the number of channels
n_channels = len(files)  # Replace with the actual number of channels
# Reshape the data into (n_channels, n_times)
min_length = min(arr.shape[0] for arr in files)

# Crop each array to the minimum length
files_target_cropped = [arr[:min_length] for arr in files]

data_reshaped = np.stack(files_target_cropped)*5e-8

# Create channel names and specify their types
ch_names = [f"{desc[0]}{desc[1]}" for desc in desc_]
ch_types = ['eeg'] * n_channels  # Adjust types as needed (e.g., 'eeg', 'ecog', 'stim')

# Create an `info` object
info = mne.create_info(ch_names=ch_names, sfreq=1375, ch_types=ch_types)

# Create an MNE RawArray object
data = mne.io.RawArray(data_reshaped, info=info)

data.plot()

data.crop(tmin=0.2)

epochs = mne.make_fixed_length_epochs(data, duration=7, preload=False)

freqs = np.arange(3, 45, 1)
n_cycles = freqs 

#welch comme script du dessus?
power_sujet = mne.time_frequency.tfr_morlet(epochs,freqs=freqs,n_cycles=n_cycles,return_itc=False,n_jobs = 13,average=True)

power_sujet.apply_baseline(mode="zscore",baseline=(None,6))
power_sujet.plot(tmin=0.2,vmin= -1,vmax = 1)