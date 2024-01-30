MESA Inlist and run_stars_extras file to reproduce the results from Sabhahit et. al (2023) - Very Massive Stars and Pair-Instability Supernovae: Mass-loss Framework for low Metallicity (https://arxiv.org/abs/2306.11785)

Mass loss for modelling the evolution of classical massive stars (10-50 Msun) is typically taken from Vink et al. (2000, 2001). However a steeper scaling of mass loss in proximity to the Eddington limit has been found in Vink et al. (2011). This study implements mass loss scaling from Vink et al. (2011) for very massive stars. In the absence of a large sample of VMS in low metallicity environments, we make use of both the concepts of transition mass loss point from Vink et al. (2012) and hydro-dynamically consistent PoWR atmosphere models to construct the new implementation. The applicable Z range for the VMS mass loss framework in this study is 0.004 to 0.0002 (SMC-like to a hundredth solar). A summary of mass loss used is
1. O stars : Vink et al. (2000, 2001)
2. VMS : Scaling from Vink et al. (2011) and absolute rate fixed at the transition mass loss point fixed from both hydro-dynamically consistent PoWR atmosphere models and Vink et al. (2012)
3. WR stars : Sander & Vink (2020)
4. Cool supergiants : de Jager (1988)

----- For high-metallicity VMS recipe -----

!! This recipe is only applicable for SMC-like or below. For high metallicity case, see 

----- Relevant files and what they do -----
MESA version r12115
1. run_iteration.py: python script to run multiple models back to back. 
2. inlist_project_H and inlist_project_H_LOGS: for core H burning
3. inlist_project_He and inlist_project_He_LOGS: for core He burning
4. /src/run_stars_extras.f: run_stars file with the low Z VMS mass loss framework

