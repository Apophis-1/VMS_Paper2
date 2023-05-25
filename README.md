# MESA Inlist and run_stars_extras file to reproduce the results from Sabhahit et. al (2023) - Very Massive Stars and Pair-Instability Supernovae: Mass-loss Framework for low Metallicity (link)
# MESA version r12115
# Applicable Z range for VMS mass loss framework - 0.008 to 0.0002. For higher Z, see VMS_paper1 repository. 

----- Relevant files and what they do -----
1. run_iteration.py: python script to run multiple models back to back. 
2. inlist_project_H and inlist_project_H_LOGS: for core H burning
3. inlist_project_He and inlist_project_He_LOGS: for core He burning
4. /src/run_stars_extras.f: run_stars file with the low Z VMS mass loss framework

