#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 17:03:55 2023

@author: gauthamsabhahit
"""

import os
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


c = 2.99792458 * (10**8)
msun = 1.99*(10**30)
lsun = 3.846 * 10**26
k = (c*msun*(10**5))/(lsun*3600*24*365)
n = 5
sigma = 5.67*10**-8
rsun = 6.955 * 10**10
G = 6.67*(10**-11)
const_k = (c*msun)/(lsun*3600*24*365)

homeDir = os.getcwd()

Z_init = 0.008                 # give the initial Z here. Values between 0.008 and 0.0002
Z_text = str(Z_init/0.02)


import fileinput
for i in [100, 200, 300]:

    j = str(i)
    for line in fileinput.FileInput("inlist_project_H_LOGS",inplace=1):
        sline=line.strip().split("=")
        if sline[0].startswith("initial_mass"):
            sline[1]= j
        elif sline[0].startswith("initial_z"):
            sline[1]= str(Z_init)
        elif sline[0].startswith("Zbase"):
            sline[1]= str(Z_init)
        elif sline[0].startswith("log_directory"):
            logs = '/Users/gauthamsabhahit/Documents/mesa-r12115/VMS_paper2/MODELS/LOGS_'+j+'_'+Z_text+'Z_sol'
            sline[1]="'"+logs+"'"
        line='='.join(sline)    
        print(line)
        
    for line in fileinput.FileInput("inlist_project_He_LOGS",inplace=1):
        sline=line.strip().split("=")
        if sline[0].startswith("Zbase"):
            sline[1]= str(Z_init)
        elif sline[0].startswith("log_directory"):
            logs = '/Users/gauthamsabhahit/Documents/mesa-r12115/VMS_paper2/MODELS/LOGS_'+j+'_'+Z_text+'Z_sol'
            sline[1]="'"+logs+"'"
        line='='.join(sline)    
        print(line)

    os.system('./clean')
    os.system('./mk')
    os.system('./rn')
