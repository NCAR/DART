Netcdf Inflation Files
======================

The filter_nml now read restart and inflation files directly from NetCDF files

Netcdf inflation files are no longer special files. DART format inflation files were always 2 copies in one file (mean
and standard devation). Taking away this special status of inflation files has the advantage that all copies (restarts,
ensemble mean, ensemble standard deviation, inflation mean, inflation sd, etc.) can all be treated the same for IO
purposes. Since there are two inflation files when reading/writing netcdf the filenames are different to DART format
restart files.

The names of the netcdf inflation files are now fixed.

**Input inflation file names**

The filter_nml option:

``inf_in_file_name = prior_inflation_ics, post_inflation_ics``

has been **deprecated** and for 1 domain filter is expecting to read:

| input_{priorinf,postinf}_mean.nc
| input_{priorinf,postinf}_sd.nc

For multiple domains filter is expecting to read:

| input_{priorinf,postinf}_mean_d01.nc
| input_{priorinf,postinf}_sd_d01.nc
| input_{priorinf,postinf}_mean_d02.nc
| input_{priorinf,postinf}_sd_d02.nc

where d0\* is the domain number.

**Output inflation file names**

The filter_nml option:

``inf_out_file_name = prior_inflation_restart, post_inflation_restart``

has been **deprecated** and for 1 domain filter is expecting to read:

| output_{priorinf,postinf}_mean.nc
| output_{priorinf,postinf}_sd.nc

For multiple domains filter is expecting to write:

| prior_inflation_restart_mean_d01
| prior_inflation_restart_sd_d01
| prior_inflation_restart_mean_d02
| prior_inflation_restart_sd_d02

where d0\* is the domain number.
