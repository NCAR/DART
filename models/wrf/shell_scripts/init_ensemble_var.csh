#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

########################################################################
#
#     init_ensemble_var.csh - script that creates perturbed initial
#                           conditions from the WRF-VAR system.
#
#     created Nov. 2007, Ryan Torn NCAR/MMM
#     modified by G. Romine and N. Collins April 2011 for use with generic 
#     wrf-dart experiments. 
#     
#     This script will prepare a set of perturbed IC files for use
#     in a wrfdart experiment. Var perturbations are added at the 
#     initial time, followed by a model advance. So, you should have
#     the initial time one assimilation window prior to the first 
#     set of observations. 
#
#
#     To run this script, the following is assumed:
#     You have already run real.exe, generating the IC/BC for the
#     initial and forecast time. These should look like:
#        - wrfinput_d01_IGDAY_IGSEC_mean
#        - wrfinput_d01_FGDAY_FGSEC_mean
#        - wrfbdy_d01_FGDAY_FGSEC_mean
#
#     these should be placed in a directory called $OUTPUT_DIR/$initial_date/
#
#     The IGDAY and IGSEC is the Gregorian day and seconds for the initial
#     time, while FGDAY and FGSEC are the Gregorian day and seconds for
#     first assimilation time (should be one assimilation window apart)
#
#     You must have da_wrfvar.exe built (assumed mpi), along with the 
#     support files of be.dat and (supplemental) bc_pert_scale
#     the executable is assumed to reside in $RUN_DIR/WRF_RUN
#
#     Further, you should have a full dart build of executables in the 
#     ${DART_DIR}/models/wrf/work directory
#
#     You must provide a script to execute a model advance starting
#     from a dart restart file
#
#     You will need to make numerous edits below to set up things to run
#     jobs on your local system. An example is provided (commented out) to
#     run on NCAR's bluefire system environment.  Look for the string: EDIT
#
########################################################################

set paramfile = init_param.csh
source $paramfile
#EDIT#
set initial_date = 2010060500
# initial_date is of the form yyyymmddhh

cd $RUN_DIR

${COPY} ${TEMPLATE_DIR}/input.nml input.nml
set gdate  = (`echo $initial_date 0h -g | ${DART_DIR}/models/wrf/work/advance_time`)
set gdatef = (`echo $initial_date ${ASSIM_INT_HOURS}h -g | ${DART_DIR}/models/wrf/work/advance_time`)
set wdate  =  `echo $initial_date 0h -w | ${DART_DIR}/models/wrf/work/advance_time`
set yyyy   = `echo $initial_date | cut -b1-4`
set mm     = `echo $initial_date | cut -b5-6`
set dd     = `echo $initial_date | cut -b7-8`
set hh     = `echo $initial_date | cut -b9-10`

${REMOVE} ${RUN_DIR}/WRF
${LINK} ${OUTPUT_DIR}/${initial_date} WRF

mkdir -p ${OUTPUT_DIR}/${initial_date}/DART
${LINK} ${OUTPUT_DIR}/${initial_date}/wrfinput_d01_${gdate[1]}_${gdate[2]}_mean wrfinput_d01

# these are the initial inflation mean and standard deviation for inflation
echo 1.0 0.6 | ${RUN_DIR}/fill_inflation_restart
${MOVE}  inflate_ics ${OUTPUT_DIR}/${initial_date}/DART/prior_inf_ic

set n = 1
while ( $n <= $NUM_ENS )

   echo "  STARTING ENSEMBLE MEMBER $n"

   set ensstring = `printf "%04d" $n`
   mkdir -p ${RUN_DIR}/mem${n}

   ${LINK} ${RUN_DIR}/WRF_RUN/* ${RUN_DIR}/mem${n}/.
   ${LINK} ${TEMPLATE_DIR}/input.nml ${RUN_DIR}/mem${n}/input.nml

   ${REMOVE} script.sed
   @ seed_array2 = $n * 1000
   cat >! script.sed << EOF
   /run_days/c\
   run_days                   = 0,
   /run_hours/c\
   run_hours                  = 0,
   /run_minutes/c\
   run_minutes                = 0,
   /run_seconds/c\
   run_seconds                = 0,
   /start_year/c\
   start_year                 = 1*${yyyy},
   /start_month/c\
   start_month                = 1*${mm},
   /start_day/c\
   start_day                  = 1*${dd},
   /start_hour/c\
   start_hour                 = 1*${hh},
   /start_minute/c\
   start_minute               = 1*00,
   /start_second/c\
   start_second               = 1*00,
   /end_year/c\
   end_year                   = 1*${yyyy},
   /end_month/c\
   end_month                  = 1*${mm},
   /end_day/c\
   end_day                    = 1*${dd},
   /end_hour/c\
   end_hour                   = 1*${hh},
   /end_minute/c\
   end_minute                 = 1*00,
   /end_second/c\
   end_second                 = 1*00,
   /analysis_date/c\
   analysis_date = \'${wdate}.0000\',
   s/PERT_SCALING/${IC_PERT_SCALE}/
   /seed_array1/c\
   seed_array1 = ${initial_date},
   /seed_array2/c\
   seed_array2 = $seed_array2 /
EOF

   sed -f script.sed ${TEMPLATE_DIR}/namelist.input >! ${RUN_DIR}/mem${n}/namelist.input
   ${LINK} ${OUTPUT_DIR}/${initial_date}/wrfinput_d01_${gdate[1]}_${gdate[2]}_mean ${RUN_DIR}/mem${n}/fg

   if ( -e ${RUN_DIR}/assim_init_${n}.csh ) then
     ${REMOVE} ${RUN_DIR}/assim_init_${n}.csh
   endif
   touch ${RUN_DIR}/assim_init_${n}.csh

   cat >> ${RUN_DIR}/assim_init_${n}.csh << EOF
#!/bin/csh

#EDIT#
########## Replace with appropriate job defs for your system, or delete

#==================================================================
#BSUB -J assim_init_${n}
#BSUB -o assim_init_${n}.%J.log
#BSUB -P ${NCAR_GAU_ACCOUNT}
#BSUB -W ${NCAR_RUNTIME}
#BSUB -q ${NCAR_QUEUE}
#BSUB -n ${NCAR_CORES}
#BSUB -x
#BSUB -R "span[ptile=${NCAR_PTILE}]"
#==================================================================

########### END of NCAR job block definitions
   
   cd ${RUN_DIR}/mem${n}

   ${MPI_EXEC} ./da_wrfvar.exe >& output.wrfvar

   ${MOVE} wrfvar_output wrfinput_d01

   ${RUN_DIR}/wrf_to_dart >& output.wrf_to_dart
   ${MOVE} dart_wrf_vector ${RUN_DIR}/filter_ic_new.${ensstring}

   ${REMOVE} wrfinput_d01

   cd $RUN_DIR
   advance_mem_restart.csh $initial_date $n $paramfile

   if ( -e ${RUN_DIR}/assim_model_state_ud.${ensstring} ) then

      ${MOVE} assim_model_state_ud.${ensstring}             ${OUTPUT_DIR}/${initial_date}/DART/filter_ic.${ensstring}
      ${MOVE} WRFOUT/wrf.out_${gdatef[1]}_${gdatef[2]}_${n} ${OUTPUT_DIR}/${initial_date}/.
      ${REMOVE} filter_ic_new.${ensstring} start_member_${n} ${RUN_DIR}/mem${n} assim_init_${n}.csh
   
   endif
EOF

# replace with an appropriate command to launch jobs on your system
     echo `eval "$JOB_SUBMIT ${RUN_DIR}/assim_init_${n}.csh"`

   @ n++

end


exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

