# aprp
A python script and some sample data are provided that demonstrate how to use the approximate partial radiative perturbation (APRP) technique to compute TOA SW radiation anomalies due to changes in individual components of the climate system. 

Sample data are from a short (2-year) period of MPI-ESM-LR using the difference between sstClimAerosol and sstClim runs.

One should difference longer periods for more robust results -- these are just for demonstrative purposes

Reference
----------
Taylor, K. E. et al. (2007), Estimating shortwave radiative forcing and response in 
    climate models, J. Clim., 20(11), 2530-2543, doi:10.1175/JCLI4143.1.  

Input
----------  
The code makes use of the following data (all provided):

| Frequency | Name | Description | Unit | File Format |
|:----------|:-----------------------------|:-------------|:------|:------------|
| monthly mean | clt | total cloud fraction | % | nc |
| monthly mean | rsdt | downwelling SW flux at the TOA | W/m^2 | nc |
| monthly mean | rsut | upwelling SW flux at the TOA | W/m^2 | nc |
| monthly mean | rsutcs | upwelling SW flux at the TOA under clear skies | W/m^2 | nc |
| monthly mean | rsuscs | upwelling SW flux at the surface under clear skies | W/m^2 | nc |
| monthly mean | rsdscs | downwelling SW flux at the surface under clear skies | W/m^2 | nc |
| monthly mean | rsds | downwelling SW flux at the surface | W/m^2 | nc |
| monthly mean | rsus | upwelling SW flux at the surface | W/m^2 | nc |

Output
----------
TOA SW flux anomalies attributable to changes in:
-surface albedo (for all-, clear-, and overcast-sky conditions)
-clouds (total change and contributions from changing cloud cover, scattering, and absorption)
-non-cloud atmosphere (e.g., from changes in water vapor, aerosols, ozone)

Each output field is size (MO,LAT,LON)


Figure Generation
----------
For the provided sample imput data, two figures are generated. These are provided in the /images/ directory
