Observation Space Diagnostics
=============================

In real applications, we don't know the truth.

We only have observations available to validate our assimilations.

DART provides a variety of plotting tools to diagnose performance using observations. 

While we know the truth in the Lorenz-96 model, it is useful to study the use of observation 
space diagnostics.

All information about observations is recorded by filter in a file using the DART observation 
sequence file format.

These files can be extended to allow complex and powerful metadata.

Before doing observation space diagnostics, the observation sequence file must be converted 
to a netcdf format with the program obs_diag.

Run obs_diag to generate a NetCDF file from the most recent assimilation.
