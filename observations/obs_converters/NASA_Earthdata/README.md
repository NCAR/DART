![DARTlogo](https://github.com/NCAR/DART/blob/Manhattan/docs/images/Dartboard7.png)

This is a quick description of the converters and utilities in this directory.

The [NASA Earthdata portal](https://earthdata.nasa.gov/)
is a front end for many (Distributed Active Archive Center) DAAC portals.

[Goddard Earth Sciences Data and Information Services Center](https://disc.gsfc.nasa.gov)

### `SMAP_L2_to_obs.f90`


### `LPRM_L3_to_obs.f90`, `AMSR_E_L2_to_obs.f90`

[LPRM_AMSRE_SOILM2: AMSR-E/Aqua surface soil moisture (LPRM) L2B V002](https://disc.gsfc.nasa.gov/datasets/LPRM_AMSRE_SOILM2_002/summary)

> AMSR-E/Aqua surface soil moisture (LPRM) L2B V002 is a Level 2 (swath) data set. 
> Its land surface parameters, surface soil moisture, land surface (skin) temperature, 
> and vegetation water content are derived from passive microwave remote sensing data from 
> the Advanced Microwave Scanning Radiometer-Earth Observing System (AMSR-E), using the 
> Land Parameter Retrieval Model (LPRM). Each swath is packaged with associated geolocation fields. 
> The data set covers the period from June 2002 to October 2011 (when the AMSR-E on the NASA EOS 
> Aqua satellite stopped producing data due to a problem with the rotation of its antenna), at the 
> spatial resolution (nominally 56 and 38 km, respectively) of AMSR-E's C and X bands (6.9 and 10.7 GHz, respectively).

### `netCDF_to_obs.f90`

The 
[Global Monthly Mean Leaf Area Index Climatology, 1981-2015](https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1653)
data may be converted with the `netCDF_to_obs` program.   Since these are monthly means,
each timestep is read and output as their own observation sequence file that has the 
date and time appended to the filename. 


