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

os.chdir("../")
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
