# VMS Mass-Loss Implementation (Low Metallicity)

MESA Inlist and run_stars_extras file to reproduce the results from [Sabhahit et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.1529S/abstract) - Very Massive Stars and Pair-Instability Supernovae: Mass-loss Framework for low Metallicity

**Links:**
- https://arxiv.org/pdf/2306.11785
- https://academic.oup.com/mnras/article/524/1/1529/7209169?searchresult=1
- https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.1529S/abstract

# Overview

Mass loss for classical massive stars (10-50 M☉) with optically-thin winds is typically taken from Vink et al. (2000, 2001), which show a shallow scaling with the Eddington parameter. However, a steeper scaling of mass loss in proximity to the Eddington limit has been found in [Vink et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011A%26A...531A.132V/abstract). This study implements a switch in the mass loss from a shallow to steeper scaling at the transition mass loss point where the single scattering limit is approximately breached as proposed in [Vink & Gräfener (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...751L..34V/abstract).

We use two different ways to estimate the transition mass loss point based on the presence or absence of observed VMSs:

**For high metallicity (GAL and LMC-like):** The transition mass-loss rate is obtained using the observed luminosities of the Of/WNh stars in the Arches cluster in our own galaxy and the 30 Dor in the LMC. For the high Z implementation, see the [VMS_paper1 repository](https://github.com/Apophis-1/VMS_Paper1) and [Sabhahit et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.3736S/abstract) - Mass-loss implementation and temperature evolution of very massive stars.

**For low metallicity (SMC-like or below):** We do not have observed individual VMS in low metallicity environments to get the transition point. In the absence of such a large sample of VMS, we make use of both the concepts of transition mass loss point from [Vink & Gräfener (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...751L..34V/abstract) and hydro-dynamically consistent PoWR atmosphere models to construct the new implementation.

**Applicable metallicity range:** Z = 0.0002 to 0.004 (a hundredth solar to SMC-like, 0.01-0.2 Z☉)

### Summary of Mass Loss Used

1. **O stars:** Vink et al. (2000, 2001)
2. **VMS:** Scaling from [Vink et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011A%26A...531A.132V/abstract) and absolute rate fixed at the transition mass loss point from both hydro-dynamically consistent PoWR atmosphere models and [Vink & Gräfener (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...751L..34V/abstract)
3. **WR stars:** Sander & Vink (2020)
4. **Cool supergiants:** de Jager (1988)

### Important Note

There are some small differences between the run_stars uploaded here on Github and the one used in the paper:

1. **Vink (2017) for stripped low mass helium stars:** This is now called for mass loss rates of stripped low mass helium stars. This is done by choosing the maximum between Vink (2017) rates and Sander & Vink (2020) rates. This is called in LINE 255.

2. **Maeder and Meynet (2000) rotation boost:** Set to DEFAULT ON in the run_stars file here. Since the mass loss boost is default on, set `mdot_omega_power = 0` in the inlist to prevent double boosting. If one wants the Friend & Abbott boost, comment out lines 197 and 285, and uncomment lines 198 and 286. Then turn on the MESA boost line.

## Relevant Files and What They Do

**MESA version:** r12115

1. **`run_iteration.py`**: Python script to run multiple models back to back
2. **`inlist_project_H`** and **`inlist_project_H_LOGS`**: MESA inlists for core H burning evolution
3. **`inlist_project_He`** and **`inlist_project_He_LOGS`**: MESA inlists for core He burning evolution
4. **`/src/run_star_extras.f`**: run_stars file with the low metallicity VMS mass loss framework implementation

## Usage

Compile the provided `run_star_extras.f` file. Modify `inlist_project_H` for core hydrogen burning or `inlist_project_He` for core helium burning to set initial mass and metallicity parameters. Run single models directly with MESA or use `run_iteration.py` to run multiple models sequentially.

**Grid parameter range:**
- **Mass range:** 10-500 M☉
- **Metallicity range:** Z = 0.0002 to 0.004 (0.01-0.2 Z☉, a hundredth solar to SMC-like)
- **Evolutionary phases:** Core hydrogen burning and core helium burning
  
## Citation

If you use this mass loss implementation in your research, please cite:
```bibtex
@ARTICLE{2023MNRAS.524.1529S,
       author = {{Sabhahit}, Gautham N. and {Vink}, Jorick S. and {Sander}, Andreas A.~C. and {Higgins}, Erin R.},
        title = "{Very massive stars and pair-instability supernovae: mass-loss framework for low metallicity}",
      journal = {\mnras},
     keywords = {stars: evolution, stars: massive, stars: mass-loss, stars: winds, outflows, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Astrophysics of Galaxies, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2023,
        month = sep,
       volume = {524},
       number = {1},
        pages = {1529-1546},
          doi = {10.1093/mnras/stad1888},
archivePrefix = {arXiv},
       eprint = {2306.11785},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.1529S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

