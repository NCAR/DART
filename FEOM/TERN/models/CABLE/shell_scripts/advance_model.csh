#!/bin/tcsh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance 
#    and copies the necessary files into the temporary directory
# 2) converts the DART output to input expected by the CABLE model
# 3) runs CABLE 
# 4) converts the CABLE model output to input expected by DART
#
# The error code from the script reflects which block it failed.
#
# Arguments are the 
# 1) process number of caller, 
# 2) the number of state copies belonging to that process, and 
# 3) the name of the filter_control_file for that process

set process = $1
set num_states = $2
set control_file = $3

echo "process      is $process"
echo "num_states   is $num_states"
echo "control_file is $control_file"

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to run CABLE
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}
echo "temp_dir is ${temp_dir}"

# Create a clean temporary directory and go there
\rm -rf  ${temp_dir}
mkdir -p ${temp_dir}
cd       ${temp_dir}

# Get the 'changing' namelist files from CENTRALDIR
# Only the namelists in CENTRALDIR have the updated information about
# the state of the model partway through an assimilation experiment.

cp -pv ../cable.nml . || exit 1
cp -pv ../input.nml . || exit 1

# FIXME ... rip out all comment lines and everything past the bang
#
# Try to ensure that the input.nml has the required value for
# dart_to_cable_nml:advance_time_present for this context.

echo '1'                       >! ex_commands
echo '/dart_to_cable_nml'      >> ex_commands
echo '/advance_time_present'   >> ex_commands
echo ':s/\.false\./\.true\./'  >> ex_commands
echo ':wq'                     >> ex_commands
( ex input.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

echo '1'                       >! ex_commands
echo '/output%averaging'       >> ex_commands
echo ':s/monthly/daily/'       >> ex_commands
echo ':wq'                     >> ex_commands
( ex cable.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

# Loop through each state
set           state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3
while($state_copy <= $num_states)
   
   set ensemble_member = `head -n $ensemble_member_line ../${control_file} | tail -n 1`
   set input_file      = `head -n $input_file_line      ../${control_file} | tail -n 1`
   set output_file     = `head -n $output_file_line     ../${control_file} | tail -n 1`

   echo "advancing ensemble member $ensemble_member ..."

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by CABLE.
   # We are going to take a CABLE netCDF restart file and simply overwrite the
   # appropriate variables. The DART output file also has the 'advance_to'
   # time - which must be communicated to the model ...
   #----------------------------------------------------------------------

   # The EXPECTED DART filter restart_out_file_name is 'dart_restart'
   # The dart_to_cable_nml:advance_time_present = .TRUE. must be set;
   # this is the purpose of the ex_commands earlier.
   # specified in cable.nml$filename%restart_in (should be a local filename)
   # gridinfo file is specified in cable.nml$filename%type (can be anything)

   set STRING=`grep "filename%restart_in" cable.nml | sed -e "s#[='\\!]# #g"`
   set CABLEIN=$STRING[$#STRING]
   set STRING=`grep "filename%restart_out" cable.nml | sed -e "s#[='\\!]# #g"`
   set CABLEOUT=$STRING[$#STRING]
   set STRING=`grep "filename%type" cable.nml | sed -e "s#[='\\!]# #g"`
   set GRIDINFO=$STRING[$#STRING]
   set STRING=`grep "casafile%cnpipool" cable.nml | sed -e "s#[='\\!]# #g"`
   set CNIPOOLFILE=$STRING[$#STRING]
   set STRING=`grep "casafile%cnpepool" cable.nml | sed -e "s#[='\\!]# #g"`
   set CNEPOOLFILE=$STRING[$#STRING]
   set STRING=`grep "gswpfile%rainf" cable.nml | sed -e "s#[='\\!]# #g"`
   set TIME_FILE=$STRING[$#STRING]

   set CABLEFILE = `printf CABLE_restart.%04d.nc     ${ensemble_member}`
   set  CNPIFILE = `printf CABLE_poolcnp_in.%04d.csv ${ensemble_member}`

   ln -sfv    ${TIME_FILE}   CABLE_time_file.nc   || exit 2
   ln -sfv    ${GRIDINFO}    CABLE_gridinfo.nc    || exit 2
   cp -fv  ../${CABLEFILE}   CABLE_restart.nc     || exit 2
   cp -fv  ../${CNPIFILE}    CABLE_poolcnp_in.csv || exit 2
   cp -fv  ../${input_file}  dart_restart         || exit 2

   ../dart_to_cable || exit 2

   #----------------------------------------------------------------------
   # Block 3: Run CABLE
   # casafile%cnpipool    = 'output/cnppool1979.csv'
   # casafile%cnpepool    = 'output/cnppool1979.csv'
   # filename%restart_in  = './restart_in.nc'
   # filename%restart_out = './restart_out.nc'
   #
   # restart_in.nc and ${CABLEIN} may be the same name ... 
   #----------------------------------------------------------------------

   mkdir -p output
   mv -fv CABLE_restart.nc     ${CABLEIN}        || exit 3
   mv -fv CABLE_poolcnp_in.csv ${CNIPOOLFILE}    || exit 3

   echo "mpi will run with the following command <${MPICMD}>"

   ${MPICMD} ../cable-mpi

   set cablestatus = $status
   if ( $cablestatus != 0 ) then
      echo "ERROR - ensemble member $ensemble_member did not complete successfully" 
      echo "ERROR - ensemble member $ensemble_member did not complete successfully" 
      exit 3 
   endif

   echo "model advanced to "
   ncdump -v time ${CABLEOUT}
   echo

   #----------------------------------------------------------------------
   # Block 4: Convert the CABLE model output to form needed by DART
   #----------------------------------------------------------------------

   # cable_to_dart reads the restart file after the model advance and writes
   # out an updated DART 'initial conditions' file. This initial conditions
   # file contains a header with the valid time of the ensuing model state.
   # The cable restart files contain the valid time of the model state.

   ln -sf ${CABLEOUT} CABLE_restart.nc

   ../cable_to_dart || exit 4

   # The (new,updated) DART restart file name is called 'dart_ics'
   # Move the updated files back to 'centraldir'
   mv -v dart_ics       ../${output_file} || exit 4
   mv -v ${CABLEOUT}    ../${CABLEFILE}   || exit 4
   mv -v ${CNEPOOLFILE} ../${CNIPOOLFILE}

   # bookkeeping

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3

end

# Change back to original directory and get rid of temporary directory
cd ..
# \rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
#\rm -rf $control_file
\mv $control_file old_control_file.$$

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

