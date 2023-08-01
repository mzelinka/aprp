# aprp

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7809085.svg)](https://doi.org/10.5281/zenodo.7809085)

A jupyter notebook is provided that demonstrates how to use the approximate partial radiative perturbation (APRP) technique to compute TOA SW radiation anomalies due to changes in individual components of the climate system for an example CMIP6 model. This notebook reads in data from two RFMIP experiments, calls the aprp.py function to perform the APRP calculations, and plots up the results. For comparison, it also performs the calculation using the formulation of Smith et al (2020), which was shown to be erroneous in Zelinka et al (2023).

APRP Output
----------
TOA SW flux anomalies attributable to changes in:

1. surface albedo (for all-, clear-, and overcast-sky conditions)
2. clouds (total change and contributions from changing cloud cover, scattering, and absorption)
3. non-cloud atmosphere (e.g., from changes in water vapor, aerosols, ozone)


References
----------
- Taylor, K. E. et al. (2007), Estimating shortwave radiative forcing and response in 
    climate models, J. Clim., 20(11), 2530-2543, [doi:10.1175/JCLI4143.1.](https://journals.ametsoc.org/doi/10.1175/JCLI4143.1)  
- Zelinka, M. D., T. Andrews, P. M. Forster, and K. E. Taylor, 2014: [Quantifying Components of Aerosol-Cloud-Radiation Interactions in Climate Models](http://onlinelibrary.wiley.com/doi/10.1002/2014JD021710/abstract), _J. Geophys. Res._, 119, 7599-7615, doi:10.1002/2014JD021710.
- Smith, C. J., et al., 2020: Effective radiative forcing and adjustments in CMIP6 models, _Atmos. Chem. Phys._, 20, 9591–9618, [doi:10.5194/acp-20-9591-2020](https://doi.org/10.5194/acp-20-9591-2020).
- Zelinka, M. D., C. J. Smith, Y. Qin, and K. E. Taylor, 2023: Comparison of Methods to Estimate Aerosol Effective Radiative Forcings in CMIP Models, _Atmos. Chem. Phys._, in press.

