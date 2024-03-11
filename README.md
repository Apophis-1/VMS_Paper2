MESA Inlist and run_stars_extras file to reproduce the results from Sabhahit et. al (2023) - Very Massive Stars and Pair-Instability Supernovae: Mass-loss Framework for low Metallicity (https://arxiv.org/abs/2306.11785)

Mass loss from optically-thin winds of classical massive stars (10-50 Msun) is typically taken from Vink et al. (2000, 2001). The mass loss from Vink et al. (2000, 2001) has a shallow dependence on the Eddington parameter. However a steeper scaling of mass loss in proximity to the Eddington limit has been found in Vink et al. (2011). This study implements mass loss scaling from Vink et al. (2011) for very massive stars. The switch from a shallow to steeper scaling occurs at the transition mass loss point as proposed in Vink et al. (2012).

We use two different ways to estimate the transition mass loss point based on the presence or absecne of observed VMS.

For high metallicity (GAL and LMC-like), the transition point is calibrated using observed Of/WNh stars in the Arches cluster in our own galaxy and the 30 Dor in the LMC. For the high Z implementation, see the VMS_paper1 repository and Sabhahit et. al (2022) - Mass-loss implementation and temperature evolution of very massive stars.

However, we do not have observed individual VMS in low metallicity environments (SMC-like or below) to calibrate the transition point. In the absence of such a large sample of VMS, we make use of both the concepts of transition mass loss point from Vink et al. (2012) and hydro-dynamically consistent PoWR atmosphere models to construct the new implementation. 

The applicable Z range for the VMS mass loss framework in this study is 0.004 to 0.0002 (SMC-like to a hundredth solar). A summary of mass loss used is
1. O stars : Vink et al. (2000, 2001)
2. VMS : Scaling from Vink et al. (2011) and absolute rate fixed at the transition mass loss point fixed from both hydro-dynamically consistent PoWR atmosphere models and Vink et al. (2012)
3. WR stars : Sander & Vink (2020)
4. Cool supergiants : de Jager (1988)

IMPORTANT NOTE: There are some small differences between the run_stars uploaded here on Github and the one used in the paper.
1. Vink (2017) is now called for mass loss rates of stripped low mass helium stars. This is done by chosing the maximum between Vink (2017) rates and Sander & Vink (2020) rates. This is called in LINE 255.
2. Maeder and Meynet (2000) rotation boost is set to DEFAULT ON in the run_stars file here. Since the mass loss boost is default on, set     "mdot_omega_power = 0"    in the inlist to prevent double boosting. If one wants the Friend & Abbott boost, comment out lines 197 and 285, and uncomment lines 198 and 286. Then turn on the MESA boost line.


----- Relevant files and what they do -----
MESA version r12115
1. run_iteration.py: python script to run multiple models back to back. 
2. inlist_project_H and inlist_project_H_LOGS: for core H burning
3. inlist_project_He and inlist_project_He_LOGS: for core He burning
4. /src/run_stars_extras.f: run_stars file with the low Z VMS mass loss framework

