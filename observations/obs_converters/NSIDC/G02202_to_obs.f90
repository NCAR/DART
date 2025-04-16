! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program G02202_to_obs

!-----------------------------------------------------------------------
!> G02202_to_obs - reads a daily netCDF file of the NOAA/NSIDC SEA ICE 
!>                 CONCENTRATION observations and writes a DART obs_seq file.
!>
!>          original - 15 April 2025 Molly Wieringa, NCAR/ASP
!>
!> "The NOAA/NSIDC Climate Data Record (CDR) of sea ice concentration 
!>  from passive microwave data is a rule-based combination of ice 
!>  concentration estimates from two well-established algorithms: the NASA 
!>  Team (NT) algorithm (Cavalieri et al. 1984) and NASA Bootstrap (BT) 
!>  algorithm (Comiso 1986). The CDR is a consistent, daily time series of 
!>  sea ice concentrations from 25 October 1978 through the most recent 
!>  processing for both the north and south polar regions. All data are on a 
!>  25 km x 25 km grid."
!>
!> Please make sure you cite the data in accordance to the agreement:
!>
!> Meier, W. N., F. Fetterer, A. K. Windnagel, J. S. Stewart, and T. Stafford. 
!> 2024. NOAA/NSIDC Climate Data Record of Passive Microwave Sea Ice Concentration, 
!> Version 5. [Indicate subset used]. Boulder, Colorado USA. NSIDC: National 
!> Snow and Ice Data Center https://doi.org/10.7265/rjzb-pf78.[Date Accessed].
!> 
!> G02202 data files are provided in netCDF format via the following link:
!> https://nsidc.org/data/g02202/versions/5#anchor-documentation. The files are 
!> housed separately for the Northern and Southern Hemispheres, but this converter
!> is written generally for either set of files. 
!> ------------------------------------------------------------------------

!  Call from required modules
! use            types_mod, only : r8

! Things to note for this converter:
!  - The G02202 data files are provided in netCDF format
!  - The files are housed separately for the Northern and Southern Hemispheres
!  - The grid is a 25x25 km polar stereographic grid, and lat/lon data is not
!    provided in the file
!    - the way YF dealt with this was to have separate files for psn25_lats/lons
!    - there are ancillary data files that have the lat/lon for each grid (NH and SH)

! Conditions to consider:
! - various masks ()

end program G02202_to_obs