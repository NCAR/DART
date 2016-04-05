#!/bin/csh
# 
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script is intentionally not 'executable' until setup_pmo.csh 
# copies it to CENTRALDIR and makes it executable.
# After everything in CENTRALDIR is confirmed to be correct, this
# script may be exectuted.
#
#-----------------------------------------------------------------------------
# These directives (for LSF) may be ignored if not running under LSF.
#
#BSUB -J jules_perfect
#BSUB -o jules_perfect.%J.log
#BSUB -P NIMG0002
#BSUB -q premium
#BSUB -n 1
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu

#-------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------

set REMOVE = 'rm -fr'

@ BAIL = 0
foreach FILE ( input.nml jules_to_dart dart_to_jules  perfect_model_obs\
               advance_model.csh jules.exe \
               urban.nml triffid_params.nml timesteps.nml prescribed_data.nml\
               pft_params.nml nveg_params.nml model_grid.nml jules_vegetation.nml\
               jules_surface_types.nml jules_surface.nml jules_soil.nml jules_snow.nml\
               jules_rivers.nml jules_radiation.nml jules_hydrology.nml\
               initial_conditions.nml imogen.nml fire.nml drive.nml crop_params.nml\
               ancillaries.nml output.nml )

   if ( ! -e $FILE ) then
      echo "$FILE is needed but not present in `pwd`"
      @ BAIL = 1
   endif

end

if ( $BAIL > 0 ) then
   echo "FATAL ERROR ... stage the missing file(s) and try again."
   echo "FATAL ERROR ... stage the missing file(s) and try again."
   exit 1
endif

./perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit -4
endif

#${MOVE} True_State.nc    ../jules_True_State.${LND_DATE_EXT}.nc
#${MOVE} obs_seq.perfect  ../jules_obs_seq.${LND_DATE_EXT}.perfect
#${MOVE} dart_log.out     ../jules_dart_log.${LND_DATE_EXT}.out

echo "`date` -- END   jules PERFECT_MODEL_OBS"

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

${REMOVE} dart_log.nml

echo "`date` -- END   GENERATE jules TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

