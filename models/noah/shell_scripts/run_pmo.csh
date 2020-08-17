#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This is an example script for how to run a perfect_model experiment.
# All the required files have been staged in CENTRALDIR and the namelists
# have been customized for the experiment. Since perfect_model_obs is a single-
# threaded program ... and we are only running one instance of NOAH ... we are
# running this on the command line.

#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
#==============================================================================

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

#==============================================================================
# Check to make sure all the required files have been staged in CENTRALDIR
#==============================================================================

set CENTRALDIR = `pwd`

@ BAIL = 0
foreach FILE ( wrfinput namelist.hrldas Noah_hrldas_beta SOILPARM.TBL \
               VEGPARM.TBL GENPARM.TBL URBPARM.TBL input.nml \
               perfect_model_obs \
               restart.nc )

   if ( ! -e $FILE ) then
      echo "$FILE is needed but not present in CENTRALDIR"
      @ BAIL = 1
   endif

end

if ( $BAIL > 0 ) then
   echo "FATAL ERROR ... stage the missing file(s) and try again."
   echo "FATAL ERROR ... stage the missing file(s) and try again."
   exit 1
endif
# run the Noah-MP 
# figure out the date of the model state, should be the same as the observation data
# run perfect_model_obs, which needs obs_seq.in, could be linked to the folder..
./perfect_model_obs

exit 0


