# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

Once MPAS is compiled for CORE=atmosphere in your MPAS_RUN directory, 
data stream files will be automatically generated as below.

streams.atmosphere
stream_list.atmosphere.output 
stream_list.atmosphere.surface 
stream_list.atmosphere.diagnostics

Thus, these default files are not provided here.

However, a sample namelist.atmosphere can be found here for cycling experiments
which are supposed to be run in a restart mode in MPAS_RUN.
For that, namelist.atmosphere should be updated as below.
(Note that these options can be turned off for the very first cycle in case 
it is run with mpas initial files rather than restart files.)
If you canot find 'config_do_DAcycling' in your default namelist.atmosphere,
you should add the parameter for your cycling runs as below.

&restart
   config_do_restart = true
   config_do_DAcycling = true
/

If config_sst_update = false in namelist.atmosphere or the input sst file does not exist, 
streams.atmosphere might be edited for input_interval="none" in the "surface" stream.

Once these files are all available and properly edited in MPAS_RUN, 
they will be used for all ensemble members throughout the cycling period, 
only being updated for config_start_time for each cycle.

At this point, the IAU option is not supported here. Please contact dart@ucar.edu if you
need help on how to use IAU in MPAS/DART.

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
