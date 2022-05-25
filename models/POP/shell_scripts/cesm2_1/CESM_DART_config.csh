#!/bin/csh
#
# Copyright 2020 University Corporation for Atmospheric Research
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License.
# Please view the License at http://www.apache.org/licenses/LICENSE-2.0
#
# ==============================================================================
#
# ---------------------
# Purpose
# ---------------------
#
# This script integrates DART with a pre-existing CESM multi-instance case.
# It must be run from a valid CASEROOT directory. If the case was created
# using one of the DART scripts, this script should be staged in the
# CASEROOT directory automatically, and DARTROOT is set at that time.
#
# POP is the only active model component.
# CESM starts and stops to allow for POP to assimilate every 24 hours.
#
# This script will build the DART executables if they are not found.
#
# There are many CESM binary files in big-endian format, and DART reads
# some of them, so you MUST compile DART accordingly e.g.,
# ifort -convert big_endian
# Contact dart@ucar.edu if you want to use another compiler.
#
# ---------------------
# How to set up the script
# ---------------------
#
# -- Ensure DARTROOT references a valid DART directory.
# -- Examine the whole script to identify things to change for your experiments.
# -- Provide any initial files needed by your run:
#       inflation
#       sampling error correction
# -- Run this script.
# -- Edit the DART input.nml that appears in the ${CASEROOT} directory.
# -- Submit the job using ${CASEROOT}/case.submit
#
# ==============================================================================
# Get the environment of the case - defines number of instances/ensemble size ...
# Each model component has their own number of instances.
# ==============================================================================

if ( ! -e ./xmlquery ) then
   echo "ERROR: $0 must be run from a CASEROOT directory".
   exit -1
endif

setenv CASE          `./xmlquery --value CASE`
setenv CASEROOT      `./xmlquery --value CASEROOT`
setenv EXEROOT       `./xmlquery --value EXEROOT`
setenv RUNDIR        `./xmlquery --value RUNDIR`
setenv NINST_OCN     `./xmlquery --value NINST_OCN`

source DART_params.csh

# ==============================================================================
# Turn on data assimilation in CESM
# 
# Data assimilation is enabled within an earth system component by changing 
# DATA_ASSIMILATION_[CPL, ATM, LND, ICE, OCN, ROF, GLC, WAV]=TRUE
# 
# ==============================================================================

./xmlchange DATA_ASSIMILATION_OCN=TRUE
./xmlchange DATA_ASSIMILATION_CYCLES=1

# ==============================================================================
# standard commands:
#
# If you are running on a machine where the standard commands are not in the
# expected location, add a case for them below.
# ==============================================================================

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

   breaksw
   default:
      # NERSC "hopper", NWSC "yellowstone", NWSC "cheyenne"
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

   breaksw
endsw

echo ""

# ==============================================================================
# make sure the required directories exist
# VAR is the shell variable name, DIR is the value
# ==============================================================================

foreach VAR ( CASEROOT DARTROOT )
   set DIR = `eval echo \${$VAR}`
   if ( ! -d $DIR ) then
      echo "ERROR: directory '$DIR' not found"
      echo "       In the setup script check the setting of: $VAR"
      exit 10
   endif
end

# ==============================================================================
# Make sure the DART executables exist or build them if we can't find them.
# The DART input.nml in the model directory IS IMPORTANT during this part
# because it defines what observation types are supported.
# ==============================================================================

foreach MODEL ( POP )
   set targetdir = ${DARTROOT}/models/${MODEL}/work
   if ( ! -x $targetdir/filter ) then
      echo ""
      echo "WARNING: executable file 'filter' not found."
      echo "         Looking for: $targetdir/filter "
      echo "         Trying to rebuild all executables for $MODEL now ..."
      (cd $targetdir; ./quickbuild.sh)
      if ( ! -x $targetdir/filter ) then
         echo "ERROR: executable file 'filter' not found."
         echo "       Unsuccessfully tried to rebuild: $targetdir/filter "
         echo "       Required DART assimilation executables are not found."
         echo "       Stopping prematurely."
         exit 20
      endif
   endif
end

# ==============================================================================
# Stage the required parts of DART in the CASEROOT directory.
# ==============================================================================

${COPY} ${DARTROOT}/models/POP/shell_scripts/${cesmtagmajor}/assimilate.csh          .
${COPY} ${DARTROOT}/models/POP/shell_scripts/${cesmtagmajor}/perfect_model.csh       .

# ==============================================================================
# Stage the DART executables in the CESM execution root directory: EXEROOT
# If you recompile the DART code (maybe to support more observation types)
# we're making a script to make it easy to install new DART executables.
# ==============================================================================

cat << EndOfText >! stage_dart_files
#!/bin/sh

# Run this script in the ${CASEROOT} directory.
# This script copies over the dart executables and POSSIBLY a namelist
# to the proper directory.  If you have to update any dart executables,
# do it in the ${DARTROOT} directory and then rerun stage_dart_files.
# If an input.nml does not exist in the ${CASEROOT} directory,
# a default one will be copied into place.
#
# This script was autogenerated by $0 using the variables set in that script.

if [[ -e input.nml ]]; then
   echo "stage_dart_files: Using existing ${CASEROOT}/input.nml"
   if [[ -e input.nml.original ]]; then
      echo "input.nml.original already exists - not making another"
   else
      ${COPY} input.nml input.nml.original
   fi

elif [[ -e ${DARTROOT}/models/POP/work/input.nml ]]; then
   ${COPY} ${DARTROOT}/models/POP/work/input.nml  input.nml
   if [[ -x update_dart_namelists ]]; then
          ./update_dart_namelists
   fi
else
   echo "ERROR: stage_dart_files could not find an input.nml.  Aborting"
   exit -99
fi

${COPY} ${DARTROOT}/models/POP/work/filter                  ${EXEROOT}
${COPY} ${DARTROOT}/models/POP/work/perfect_model_obs       ${EXEROOT}
${COPY} ${DARTROOT}/models/POP/work/fill_inflation_restart  ${EXEROOT}

exit 0

EndOfText
chmod 0755 stage_dart_files

./stage_dart_files  || exit -8

# ==============================================================================
# Ensure the DART namelists are consistent with the ensemble size,
# suggest settings for num members in the output diagnostics files, etc.
# The user is free to update these after setup and before running.
# ==============================================================================

cat << EndOfText >! update_dart_namelists
#!/bin/sh

# this script makes certain namelist settings consistent with the number
# of ensemble members built by the setup script.
# this script was autogenerated by $0
# using the variables set in that script

# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output observation sequence files contain ALL ensemble members.

ex input.nml <<ex_end
g;ens_size ;s;= .*;= ${NINST_OCN};
g;num_output_obs_members ;s;= .*;= ${NINST_OCN};
wq
ex_end

exit 0

EndOfText
chmod 0755 update_dart_namelists

./update_dart_namelists || exit 70

#=========================================================================
# Stage the files needed for SAMPLING ERROR CORRECTION - even if not
# initially requested. The file is static, small, and may be needed later.
#
# If it is requested and is not present ... it is an error.
#
# The sampling error correction is a lookup table.  A selection of common
# ensemble sizes should be found in the file named below.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

if ( $NINST_OCN == 1 ) then
    ./xmlchange DATA_ASSIMILATION_SCRIPT=${CASEROOT}/perfect_model.csh
else
    ./xmlchange DATA_ASSIMILATION_SCRIPT=${CASEROOT}/assimilate.csh
    
    set SAMP_ERR_DIR = assimilation_code/programs/gen_sampling_err_table/work
    set SAMP_ERR_FILE = ${DARTROOT}/${SAMP_ERR_DIR}/sampling_error_correction_table.nc
 
    if (  -e   ${SAMP_ERR_FILE} ) then
        ${COPY} -f ${SAMP_ERR_FILE} ${RUNDIR}  || exit 75
        if ( $NINST_OCN < 3 || $NINST_OCN > 200 ) then
            echo ""
            echo "ERROR: sampling_error_correction_table.nc handles ensemble sizes 3...200."
            echo "ERROR: Yours is ${NINST_OCN}"
            echo ""
            exit 75
        endif
    else
        set list = `grep sampling_error_correction input.nml | sed -e "s/[=\.,]//g`
        if ($list[2] == "true") then
            echo ""
            echo "ERROR: No sampling_error_correction_table.nc file found ..."
            echo "ERROR: the input.nml:assim_tool_nml:sampling_error_correction"
            echo "ERROR: is 'true' so this file must exist."
            echo ""
            exit 80
        endif
    endif
endif

# ==============================================================================
# INFLATION : Initial setup for (default) adaptive state-space prior inflation
# ==============================================================================

# The initial inflation values may be specified in two ways. One is via the
# namelist, the other is to run 'fill_inflation_restart' to create netCDF
# inflation files. Either of these methods should ONLY be run on the first cycle.
# The presence of the pop_inflation_cookie indicates the first cycle.

if ( ${NINST_OCN} > 1 ) then
   date >! ${RUNDIR}/pop_inflation_cookie
endif

# ==============================================================================
# What to do next
# ==============================================================================


cat << EndOfText >! DART_instructions.txt

-------------------------------------------------------------------------

Check the DART configuration:

1) The default behavior is to _not_ invoke DART and simply run CESM.
   We recommend that you make sure this works before proceeding.

2) To turn DART on or off, while in ${CASEROOT}  execute:
   ./xmlchange --value DATA_ASSIMILATION_OCN=[TRUE,FALSE] 

3) Modify what you need to in the DART namelist file, i.e. ${CASEROOT}/input.nml

4) If you have recompiled any part of the DART system, 'stage_dart_files'
   will copy them into the correct places.

5) If you stage your own inflation files, make sure you read the "INFLATION" section
   in ${CASEROOT}/assimilate.csh. At the very least, copy your inflation files
   in ${RUNDIR} into the appropriate names:
      ${CASE}.pop.output_priorinf_mean.nc
      ${CASE}.pop.output_priorinf_sd.nc
      ${CASE}.pop.output_postinf_mean.nc 
      ${CASE}.pop.output_postinf_sd.nc   
   and remove the ${RUNDIR}/pop_inflation_cookie file. If that file remains,
   your initial inflation values will come from the namelist! The staged files
   should be 'output' as the logic of assimilate.csh renames the output of the
   previous assimilation cycle to be reused as the input for next cycle.  

   If assimilate.csh does not find inflation files, it will call fill_inflation_restart
   to create some from the inflation values in input.nml.

6) Make sure the observation directory names in assimilate.csh or perfect_model.csh
   matches the one on your system.

7) Submit the CESM job in the normal way.

8) You can use ${CASEROOT}/stage_cesm_files to stage an ensemble of files
   to restart a run at a date for which you have a restart set.

9) You can run multiple forecast/assimilation cycles in a single job by
   setting the DATA_ASSIMILATION_CYCLES xml variable. To perform 4 cycles
   in a single job - specify:
   ./xmlchange --value DATA_ASSIMILATION_CYCLES=4
   make sure you give the job enough wallclock to finish all 4 cycles,
   and be prepared for your RUNDIR to be very full.

-------------------------------------------------------------------------

EndOfText

cat DART_instructions.txt

exit 0

