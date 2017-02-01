#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

set SNAME = $0
set clobber

if ($#argv == 1) then

      setenv LOGDIR $1
      echo "Running $SNAME:t and archiving results to $LOGDIR"
  
else
      echo " "
      echo "usage: $SNAME:t destinationdir"
      echo " "
      echo "This script runs the 'workshop_setup.csh' script for a wide range of models"
      echo "and archives the results in 'destinationdir'. If 'destinationdir' does not"
      echo "exist, it is created. If it does exist, any duplicate contents will be"
      echo "overwritten. This must be run from the top-level 'DART' directory."
      echo " "
      echo "This is a pretty verbose process, so if you are logging the output,"
      echo "make sure you have plenty of space:"
      echo " "
      echo "./$SNAME:t archdir |& tee DART_test.log"
      echo " "
      echo "can easily result in a 2 MB log file"
      echo " "
      exit 1
endif

if ( ! -d models/lorenz_96 ) then
   echo "$SNAME:t must be run from the top-level"
   echo "DART directory -- please try again."
   exit 2
else
   set DARTHOME = `pwd`
endif

setenv ARCHIVEDIR ${DARTHOME}/${LOGDIR}
if ( ! -d $ARCHIVEDIR ) then
   echo "Creating $ARCHIVEDIR"
   mkdir -pv $ARCHIVEDIR
endif

echo "The top-level DART directory (DARTHOME) is $DARTHOME"

set REMOVE = 'rm -fv'
set   MOVE = 'mv -f'

#----------------------------------------------------------------------
# Run the workshop_setup demo for a wide range of models.
#
# PBL_1d             needs -r8 compile flag
# MITgcm_annulus     needs location_mod changes and model_mod changes
#----------------------------------------------------------------------

@ makenum  = 1
@ modelnum = 101
foreach MODEL ( template ikeda 9var lorenz_63 lorenz_84 lorenz_96 \
                lorenz_04 lorenz_96_2scale forced_lorenz_96 \
                bgrid_solo pe2lyr null_model cam wrf )

   echo "----------------------------------------------------------"
   echo "Running $MODEL at "`date`
   echo ""

   cd ${DARTHOME}/models/${MODEL}/work

#  ${REMOVE}           obs_seq.in obs_seq.out obs_seq.final \
   ${REMOVE} input.nml obs_seq.in obs_seq.out obs_seq.final \
             perfect_ics filter_ics \
             True_State.nc Prior_Diag.nc Posterior_Diag.nc

   svn update

   csh workshop_setup.csh 

   set runstat = $status

   if ($runstat > 0) then
     echo "ERROR: model $MODEL workshop status was $runstat"
     echo "ERROR: not good."
     exit $modelnum
   endif

   if ( -e obs_seq.final ) then

      switch ( $MODEL )
         case ikeda:
            echo "$MODEL has identity obs - no point running obs_diag."
            breaksw
         case 9var:
            echo "$MODEL has identity obs - no point running obs_diag."
            breaksw
         case lorenz_63:
            echo "$MODEL has identity obs - no point running obs_diag."
            breaksw
         case lorenz_84:
            echo "$MODEL has identity obs - no point running obs_diag."
            breaksw
         case lorenz_04:
            echo "$MODEL has identity obs - no point running obs_diag."
            breaksw
         case lorenz_96_2scale:
            echo "$MODEL has identity obs - no point running obs_diag."
            breaksw
         case default:
            ./obs_diag
            breaksw
      endsw

   endif

   #------------------------------------------------------------------
   # Save the output to a directory for comparison
   #------------------------------------------------------------------

   mkdir -pv ${ARCHIVEDIR}/${MODEL}/work

   ${MOVE} dart_log.out        ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} True_State.nc       ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} perfect_restart     ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} obs_seq.out         ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} Prior_Diag.nc       ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} Posterior_Diag.nc   ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} obs_seq.final       ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} filter_restart      ${ARCHIVEDIR}/${MODEL}/work

   ${MOVE} *ges_times*dat      ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} *anl_times*dat      ${ARCHIVEDIR}/${MODEL}/work
   ${MOVE} ObsDiagAtts.m       ${ARCHIVEDIR}/${MODEL}/work

   #------------------------------------------------------------------
   # clean up and move on
   #------------------------------------------------------------------
   ${REMOVE} *.o *.mod Makefile .cppdefs input.nml*default 
   ${REMOVE} preprocess filter perfect_model_obs create_fixed_network_seq 
   ${REMOVE} create_obs_sequence obs_diag merge_obs_seq wakeup_filter

   @ makenum  = $makenum  + 1
   @ modelnum = $modelnum + 1

end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

