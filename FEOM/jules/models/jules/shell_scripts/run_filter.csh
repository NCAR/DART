#!/bin/csh
# 
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script is designed to be submitted as a batch job but may be run from
# the command line on some systems. Be aware that if you are
# running interactively, the manner in which you invoke filter matters.
#
#-----------------------------------------------------------------------------
#
#BSUB -J jules_filter
#BSUB -o jules_filter.%J.log
#BSUB -P P86850054
#BSUB -q premium
#BSUB -n 1
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_HOSTS) then

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv MPI_RUN_CMD mpirun.lsf

else

   #-------------------------------------------------------------------
   # You can also run this interactively
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     jules_filter
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv MPI_RUN_CMD ''

endif

#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $MYHOST"
echo "${JOBNAME} ($JOBID) started   at "`date`
echo

#-------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------
set REMOVE = 'rm -fr'

@ BAIL = 0
foreach FILE ( input.nml jules_to_dart dart_to_jules filter\
               advance_model.csh jules.exe \
               urban.nml triffid_params.nml timesteps.nml prescribed_data.nml\
               pft_params.nml nveg_params.nml model_grid.nml jules_vegetation.nml\
               jules_surface_types.nml jules_surface.nml jules_soil.nml jules_snow.nml\
               jules_rivers.nml jules_radiation.nml jules_hydrology.nml\
               initial_conditions.nml imogen.nml fire.nml drive.nml crop_params.nml\
               ancillaries.nml output.nml )

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

${MPI_RUN_CMD} ./filter

if ($status != 0) then
   echo "ERROR ... DART died in 'filter' ... ERROR"
   echo "ERROR ... DART died in 'filter' ... ERROR"
   exit -4
endif

#${MOVE} True_State.nc    ../jules_True_State.${LND_DATE_EXT}.nc
#${MOVE} obs_seq.perfect  ../jules_obs_seq.${LND_DATE_EXT}.perfect
#${MOVE} dart_log.out     ../jules_dart_log.${LND_DATE_EXT}.out

echo "`date` -- END   jules   filter"

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

${REMOVE} dart_log.nml

echo "`date` -- END  ASSIMILATE" 

exit 0

