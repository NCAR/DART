#!/bin/tcsh
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Advance a single instance of the COAMPS model by:
#  - copying the analysis from the DART directory to the run directory 
#  - converting the DART state vector to a COAMPS restart file
#  - advancing the COAMPS model
#  - converting the COAMPS restart file to a DART state vector
#  - copying the forecast out of the run directory to the DART directory 
#
# This script really shouldn't be called by itself, but rather by the
# advance_group.csh script that handles advancing a *group* of model 
# instances by repeatedly executing this script.
#
# The command line arguments specify a working directory, the number of
# the ensemble member, the number of processes to use for advancing the
# model, and the name of a named pipe that we use for communicating back
# to the script that handles integrating several model instances.
######

set PBS_O_WORKDIR = $1
set element = $2
set NPROCS = $3
set PIPENAME = $4

# Import the job-specific resource commands
source ${PBS_O_WORKDIR}/job_setup.csh

# Where the main programs live 
set COAMPS_BIN_DIR = /home/timw/COAMPS/coamps3.1.1/bin
set COAMPS_DART_DIR = /home/timw/sandbox/DART/models/coamps

# Script to do the namelist stripping
set NL_STRIP=${COAMPS_DART_DIR}/shell_scripts/strip_namelist.pl

# The COAMPS ensemble configuration has directories associated with
# the ensemble member number, thus member numbers with digits less
# than the maximum must be zero-padded.  For ease of use, just force
# this to a fixed digit length for now.
set digits      = 5
set element_num = `printf "%0${digits}d" $element`
set ensdir      = $element_num

# Where we'll do the actual model processing
set COAMPS_ENS = /net/ds-0b/export/ensembles/test/jamaica_upgrade/${ensdir}
echo ${COAMPS_ENS}
cd ${COAMPS_ENS}

# Get the date-time group
set dtg = `grep cdtg namelist | awk -F"'" '{print $2}'`

# Move into the data directory and start bringing data files in
cd ${COAMPS_ENS}/data

# Move the dart initial condition file to the COAMPS directory: do
# this as a copy since if the integration fails, we'll retry and
# want to have the ic file available.
set ic_file = `printf assim_model_state_ic.%04d ${element_num}`
cp ${PBS_O_WORKDIR}/${ic_file} ${COAMPS_ENS}/data/dart_vector

# Remove the converter-specific details from the COAMPS namelist
${NL_STRIP} ../namelist convert.vars convert

# Convert the dart initial condition file to a COAMPS restart file
rm -f dart.kstart
${PBS_O_WORKDIR}/trans_dart_to_coamps
if ($? != 0) then
    echo "Unsuccessful translation from DART of member $element!"
    echo "${element}.ERROR" > $PIPENAME
    exit 1
endif

# The DART initial condition file contains information about the current
# and target forecast times - we need to update the COAMPS namelist to
# reflect this information.  This currently works by having the translate
# program dump the information in a dart.kstart file that we read here.
@ ktaust_hour = `awk -F, '{print $1}' dart.kstart`
@ ktaust_min  = `awk -F, '{print $2}' dart.kstart`
@ ktaust_sec  = `awk -F, '{print $3}' dart.kstart`
@ ktauf_hour  = `awk -F, '{print $4}' dart.kstart`
@ ktauf_min   = `awk -F, '{print $5}' dart.kstart`
@ ktauf_sec   = `awk -F, '{print $6}' dart.kstart`

# Update the time variables in the COAMPS namelist - need to do both 
# ktau(st|f) and ksavea.  Set ksavea (the time interval to dump a restart
# file) to be equal to ktauf (the end of the forecast interval) to force a 
# restart file write only at the end of a run.
# These replaces work in two steps using placeholders - it would be nice
# to get it working where they could do it in one pass.
sed -i.bak -e 's/\(\s*ktaust\s*=\).*,\s*$/\1 NEW_KTAUSTH,NEW_KTAUSTM,NEW_KTAUSTS,/' -e "s/NEW_KTAUSTH/${ktaust_hour}/" -e "s/NEW_KTAUSTM/${ktaust_min}/" -e "s/NEW_KTAUSTS/${ktaust_sec}/" ${COAMPS_ENS}/namelist

sed -i.bak -e 's/\(\s*ktauf\s*=\).*,\s*$/\1 NEW_KTAUFH,NEW_KTAUFM,NEW_KTAUFS,/' -e "s/NEW_KTAUFH/${ktauf_hour}/"  -e "s/NEW_KTAUFM/${ktauf_min}/"  -e "s/NEW_KTAUFS/${ktauf_sec}/" ${COAMPS_ENS}/namelist
 
sed -i.bak -e 's/\(\s*ksavea\s*=\).*,\s*$/\1 NEW_KSAVEAH,NEW_KSAVEAM,NEW_KSAVEAS,/' -e "s/NEW_KSAVEAH/${ktauf_hour}/"  -e "s/NEW_KSAVEAM/${ktauf_min}/"  -e "s/NEW_KSAVEAS/${ktauf_sec}/" ${COAMPS_ENS}/namelist

# Run the COAMPS forecast
cd ${COAMPS_ENS}
echo "Running COAMPS ensemble member $element"
mpiexec -pernode -n ${NPROCS} ${COAMPS_BIN_DIR}/atmos_forecast.exe
set mpi_status = $?
if ($mpi_status != 0) then
    echo "Unsuccessful completion of member $element!"
    echo "${element}.ERROR" > $PIPENAME
    exit 1
endif

# Now that the we've run the model, move into the COAMPS data directory 
# and move things back to the DART directory. 
cd ${COAMPS_ENS}/data
rm dart_vector

# Convert the resulting COAMPS file to a dart file
${NL_STRIP} ../namelist convert.vars convert
${PBS_O_WORKDIR}/trans_coamps_to_dart
if ($? != 0) then
    echo "Unsuccessful translation to DART of member $element!"
    echo "${element}.ERROR" > $PIPENAME
    exit 1
endif

# Move the dart file back to the working directory -  only now can we
# safely delete the initial condition file.
set ud_file = `printf assim_model_state_ud.%04d ${element_num}`
mv dart_vector $PBS_O_WORKDIR/${ud_file}
rm ${PBS_O_WORKDIR}/${ic_file} 

cd ${COAMPS_ENS}

# Clean up and signal completion to the parent script
rm *.bak
echo "${element}.RAN" > $PIPENAME

exit 0


