#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#---------------------------------------------------------------------------------------------------
# Script prepare_ideal_IC.csh
# Original script by David Dowell
# Modified by Altug Aksoy, 01/06/05
#
# Purpose: Creates an initial ensemble using the WRF ideal.exe command.
#
# Needed:  1. ideal.input.filter in this directory (ideal.input will be linked to this in run directory)
#          2. if also running control, ideal.input.control in this directory
#          3. similar to 1 & 2, namelist.input.filter and namelist.input.control in this directory
#          4. Sounding files (perturbed/unperturbed) must be in the sounding directory
#          5. wrf_to_dart (executable) in this directory
#          6. all files necessary to run ideal.exe in run directory (including base/perturbed sounding)
#
# Requires 3 command-line input arguments: m -> 0: no control ics file, 1: control ics file
#                                                  Running the control option produces the input
#                                                  file "perfect_ics"; control run is based on the
#                                                  "base sounding", i.e. file "input_sounding"
#                                          n -> 0: no enmseble ics file, 1: ensemble ics file
#                                                  Running the ensemble option produces the input
#                                                  file "filter_ics"
#                                          p -> 0: no perturbed soundings, 1:perturbed soundings

#---------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------
# User-supplied script variables.
#---------------------------------------------------------------------------------------------------

set ES =  50                                                       # ensemble size

setenv RUN_DIR /users/romine/assim_dir/jan2010/rad_regression/IC   # ensemble (and control) run directory
                                                                   # (this is from where all WRF files will
                                                                   # be copied to a temporary directory)

setenv SNDG_DIR ${RUN_DIR}/sounding_perturbation                  # location of sounding files (perturbed
                                                                     # or unperturbed)

#---------------------------------------------------------------------------------------------------
# End of user specifications.
#---------------------------------------------------------------------------------------------------


# Check for arguments
# Note that for argument n=1, the additional input file perfect_ics is generated
if ($#argv != 3) then
   echo
   echo "Usage of this script is: csh prepare_ideal_IC.csh m n p"
   echo "                         m=0: skip generation of the control input (i.e., perfect_ics)"
   echo "                         m=1: generate control input (i.e., perfect_ics)"
   echo "                         n=0: skip generation of the ensemble input (i.e., filter_ics)"
   echo "                         n=1: generate ensemble input (i.e., filter_ics)"
   echo "                         p=0: no perturbed sounding"
   echo "                         p=1: perturbed sounding"
   echo
   exit 1
else
   set DC = $1
   set DE = $2
   set DP = $3
endif

if ( (! $DC ) && (! $DE ) ) then
   echo "Both control and ensemble options set to zero. Doing nothing..."
   echo
   exit 1
endif

setenv WORK_DIR `pwd`
#setenv RUN_ROOT ${RUN_DIR}/..
#setenv RUN_TEMP ${RUN_ROOT}/prepare_ideal_temp
setenv RUN_TEMP ${RUN_DIR}/prepare_ideal_temp

# Check first if wrf_to_dart exists in the current (working) directory
# It is needed in all conversions from WRF NetCDF format to DART format
if ( (! -e ../wrf_to_dart) || (! -e ./input.nml) ) then
   echo
   echo Executable "wrf_to_dart" or its required input.nml not found in ${WORK_DIR}
   echo
   exit
endif

# Also check if the temporary directory exists
# If it exists, it will be removed first
if ( -e ${RUN_TEMP} ) then
   \rm -fr ${RUN_TEMP}
endif

# Now links will be created to ${RUN_DIR} directory so first check if
# all necessary files exist in there
if ( (! -e ${RUN_DIR}/../WRF_RUN/ideal.exe) || (! -e ${RUN_DIR}/../WRF_RUN/gribmap.txt)  ) then
   echo
   echo Some of files ideal.exe, gribmap.txt
   echo do not exist in directory ${RUN_DIR}/../WRF_RUN/
   echo
   exit
endif
if (! $DP ) then
   if (! -e ${SNDG_DIR}/input_sounding) then
      echo File input_sounding does not exist in directory ${SNDG_DIR}
      echo
      exit
   endif
else
   if ( ($DC) && (! -e ${SNDG_DIR}/input_sounding) ) then
      echo File input_sounding does not exist in directory ${SNDG_DIR} for the control run.
      echo
      exit
   endif
   set NC = 1
   set file_missing = 0
   while ( $NC <= $ES )
      if (! -e ${SNDG_DIR}/input_sounding$NC) then
         set file_missing = 1
      endif
      @ NC ++
   end
   if ( $file_missing ) then
#      echo Some of perturbed sounding files (1-$ES)
      echo Some of perturbed sounding files 
      echo do not exist in directory ${SNDG_DIR}
      echo
      exit
   endif
endif

# Make links now
mkdir ${RUN_TEMP}
cd    ${RUN_TEMP}
ln -sf ${RUN_DIR}/../WRF_RUN/ideal.exe ideal.exe
ln -sf ${RUN_DIR}/../WRF_RUN/gribmap.txt gribmap.txt

#---------------------------------------------------------------------------------------------------
# Create control run (runs WRF ideal.exe, converts WRF to DART). 
#---------------------------------------------------------------------------------------------------

if ( $DC ) then
   echo
   echo Generating WRF initial conditions and DART "perfect_ics" file for the control run
   echo

   # Link to base input sounding
   ln -sf ${SNDG_DIR}/input_sounding input_sounding

   # First check if ideal.input.control & namelist.input.control exist so they can be linked to
#   if ( ! -e ${WORK_DIR}/ideal.input.control ) then
#      echo
#      echo File "ideal.input.control" not found in ${WORK_DIR}
#      echo
#      cd ${WORK_DIR}
#      exit
#   endif

   if ( ! -e ${WORK_DIR}/namelist.input.control ) then
      echo
      echo File "namelist.input.control" not found in ${WORK_DIR}
      echo
      cd ${WORK_DIR}
      exit
   endif

#   # Now link to ideal.input.control & namelist.input.control so that ideal.exe can run
#   ln -sf ${WORK_DIR}/ideal.input.control ideal.input
   ln -sf ../namelist.input.control namelist.input

   # Run ideal.exe to produce initial conditions wrfinput_d01
   \rm -f wrfinput_d01 show_domain_0000 rsl.out.0000 rsl.error.0000
   ideal.exe

   # Convert WRF initial condition file wrfinput_d01 (netcdf) into DART format
   cd ${WORK_DIR}
   \rm -f perfect_ics
   \rm -f wrfinput_d01
   \rm -f out.wrf_to_dart
   cp ${RUN_TEMP}/wrfinput_d01 wrfinput_d01
   ../wrf_to_dart >& out.wrf_to_dart
   cat dart_wrf_vector >> perfect_ics
   \rm -f dart_wrf_vector
   mv wrfinput_d01 wrfinput_d01_control

   # Go back to run directory to generate ensemble member files
   cd ${RUN_TEMP}
   \rm -f ideal.input
   \rm -f namelist.input

else
   echo
   echo Initial condition generation for the control run skipped...
   echo

endif


#---------------------------------------------------------------------------------------------------
# Create ensemble members (runs WRF ideal.exe, converts WRF to DART). 
#---------------------------------------------------------------------------------------------------

if ( $DE ) then
   echo
   echo Generating WRF initial conditions and DART "filter_ics" file for the ensemble run
   echo $ES members will be generated
   echo

   # Go to run directory
   cd ${RUN_TEMP}
   \rm -f ideal.input
   \rm -f namelist.input

#   # First check if ideal.input.filter & namelist.input.filter exist so they can be linked to
#   if ( ! -e ${WORK_DIR}/ideal.input.filter ) then
#      echo
#      echo File "ideal.input.filter" not found in ${WORK_DIR}
#      echo
#      cd ${WORK_DIR}
#      exit
#   endif

#   if ( ! -e ${WORK_DIR}/namelist.input.filter ) then
#      echo
#      echo File "namelist.input.filter" not found in ${WORK_DIR}
#      echo
#      cd ${WORK_DIR}
#      exit
#   endif

   # Now link to ideal.input.filter & namelist.input.filter so that ideal.exe can run
#   ln -sf ${WORK_DIR}/ideal.input.filter ideal.input
   ln -sf ${WORK_DIR}/../namelist.input namelist.input
   ln -sf ${WORK_DIR}/init_ideal.ncl init_ideal.ncl

\rm -f ${WORK_DIR}/filter_ics
touch ${WORK_DIR}/filter_ics

   # Loop over the ensemble members
   set NC = 1
   set runfile="nclrun.out"
   rm  ${WORK_DIR}/$runfile
   while ( $NC <= $ES )
      echo
      echo Working on ensemble member $NC
      echo

      # Link to appropriate (base or perturbed) input sounding
      if (! $DP ) then
         echo Based on base sounding data ${SNDG_DIR}/input_sounding
         ln -sf ${SNDG_DIR}/input_sounding input_sounding
      else
         echo Based on perturbed sounding data ${SNDG_DIR}/input_sounding$NC
         ln -sf ${SNDG_DIR}/input_sounding$NC input_sounding
      endif

      # Run ideal.exe to produce initial conditions wrfinput_d01
      \rm -f wrfinput_d01 show_domain_0000 rsl.out.0000 rsl.error.0000
       ideal.exe
       set cmd="/usr/local/bin/ncl 'member_num=$NC' init_ideal.ncl"
       echo $cmd
#       $cmd
       touch ${WORK_DIR}/$runfile
       echo "$cmd" >> ${WORK_DIR}/$runfile
       chmod +x ${WORK_DIR}/$runfile
       ${WORK_DIR}/$runfile
       \rm -f ${WORK_DIR}/$runfile    

      # Convert WRF member initial condition file wrfinput_d01 (netcdf) files into DART format
      cd ${WORK_DIR}
      \rm -f wrfinput_d01
      \rm -f out.wrf_to_dart
      ln -sf ${RUN_TEMP}/wrfinput_d01 wrfinput_d01
      ../wrf_to_dart >& out.wrf_to_dart
      cat dart_wrf_vector >> ${WORK_DIR}/filter_ics
      \rm -f dart_wrf_vector
      cd ${RUN_TEMP}

      @ NC ++
   end

   cd ${WORK_DIR}

else
   cd ${WORK_DIR}
   echo
   echo Initial condition generation for the ensemble run skipped...
   echo

endif


echo
echo Script prepare_ideal_IC.csh completed...
echo 

exit 0


