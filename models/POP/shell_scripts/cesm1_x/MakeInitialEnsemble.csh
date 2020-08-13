#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#=============================================================================
# So the background is that I want 40 restart files to create an initial
# conditions file with 40 ensemble members. 
#
# For lack of a better idea: I ran the model starting from 1996 01 01 ??Z 
# for 20 days with a 900 second timestep and a dumpFreq of 12 hours.
# This actually generates 41 restart files, but who's counting.
# I just deleted the LAST restart before running this script.
#
# We're just going to make the first restart into ensemble member 1. 
# Second restart will be ensemble member 2 ... I'm just not going to use
# the last restart. This may not be the greatest idea - I don't know
# how much spread we will have to start. There must be SOME ...
# depends on the nonlinearity of the model.
#
# Repeatedly call pop_to_dart to turn each restart into a 
# DART initial condition file. Each DART initial conditions file 
# starts with a timestamp that reflects the valid time for the ensuing 
# model state. To generate an ensemble of initial conditions files
# from successive restart files ... we will need to overwrite this
# timestamp so that ALL the DART initial conditions files reflect 
# the same timestamp rather than the true timestamp.
#
# The call to restart_file_tool changes the time tag in the file
# to be whatever is in the namelist. In this case, we want the time
# tag to be 1996 01 01 00Z ... in DART-speak ... 144270 days 0 seconds.
# tag to be 2000 01 01 00Z ... in DART-speak ... 145731 days 0 seconds.
#
# TJH - 28 July 2008

# WARNING: You should run this in an 'empty' directory
# WARNING: You should run this in an 'empty' directory
# WARNING: You should run this in an 'empty' directory

# Get input for this execution

set RESTARTDIR = /ptmp/thoar/POP/CCSM_x3/rest
set     POPDIR = /fs/image/home/thoar/SVN/DART/models/POP/input
set    DARTEXE = /fs/image/home/thoar/SVN/DART/models/POP/work
set    DARTNML = /fs/image/home/thoar/SVN/DART/models/POP/work/input.nml

# Check to see if all the pieces exist

set checkstat = 0

if ( -d ${RESTARTDIR} ) then
   echo "Using restart files from ${RESTARTDIR}"
else
   echo "ERROR:restart files directory ${RESTARTDIR} does not exist."
   set checkstat = 1
endif

if ( -f ${DARTNML} ) then
   echo "Using ${DARTNML}"
else
   echo "ERROR: ${DARTNML} does not exist."
   set checkstat = 1
endif

if ( -f ${DARTEXE}/pop_to_dart ) then
   echo "Using DART pop_to_dart from ${DARTEXE}"
else
   echo "ERROR: ${DARTEXE}/pop_to_dart does not exist."
   echo " Build ${DARTEXE}/pop_to_dart"
   set checkstat = 1
endif

if ( -f ${DARTEXE}/restart_file_tool ) then
   echo "Using DART restart_file_tool from ${DARTEXE}"
else
   echo "ERROR: ${DARTEXE}/restart_file_tool does not exist."
   echo " Build ${DARTEXE}/restart_file_tool"
   set checkstat = 1
endif

if ( $checkstat > 0 ) then
   echo "ERROR: required input failure."
   exit 9
endif

#-------------------------------------------------------------------------
# The whole operation relies on the fact that there is a directory
# somewhere that has the restart files and the corresponding input
# files: data, data.cal - as well as (potentially) another directory
# that has the DART namelists and executables.
#
# We will create a little subdirectory off the DARTEXE directory to
# contain the DART restart files, etc.
#-------------------------------------------------------------------------
# pop_to_dart needs files:  

cp ${DARTNML}                            .
cp ${POPDIR}/pop_in                      .
cp ${POPDIR}/horiz_grid.gx3v5.r8ieee.le  .
cp ${POPDIR}/topography.gx3v5.i4ieee.le  .
cp ${POPDIR}/vert_grid.gx3v5             .

#-------------------------------------------------------------------------
# Ensure the namelists have the right values.
#-------------------------------------------------------------------------
# Need to modify rest of input.nml for test run
# Essentially, we want the namelist block to look like:
#&restart_file_tool_nml
#   input_file_name              = "dart_input",
#   output_file_name             = "dart_output",
#   ens_size                     =  1,
#   single_restart_file_in       = .true.,
#   single_restart_file_out      = .true.,
#   write_binary_restart_files   = .true.,
#   overwrite_data_time          = .true.,
#   new_data_days                =  145731,
#   new_data_secs                =       0,
#   input_is_model_advance_file  = .false.,
#   output_is_model_advance_file = .false.,
#   overwrite_advance_time       = .false.,
#   new_advance_days             =  -1,
#   new_advance_secs             =  -1   /

echo ':0'                               >! ex_commands
echo '/restart_file_tool_nml'           >> ex_commands
echo '/write_binary_restart_files'      >> ex_commands
echo ':s/.false./.true./'               >> ex_commands
echo '/overwrite_data_time'             >> ex_commands
echo ':s/.false./.true./'               >> ex_commands
echo '/new_data_days'                   >> ex_commands
echo ':s/-1/145731/'                    >> ex_commands
echo '/new_data_secs'                   >> ex_commands
echo ':s/-1/0/'                         >> ex_commands
echo ':wq'                              >> ex_commands

( ex input.nml < ex_commands) >& /dev/null
\rm -f ex_commands

#-------------------------------------------------------------------------
# Loop over all the restart files.
#-------------------------------------------------------------------------

@ member = 0

foreach FILE ( ${RESTARTDIR}/*/cx3.dart.001.pop.r.*.nc )

   #----------------------------------------------------------------------
   # link each file to the expected name, and make sure the pop pointer
   # file has the same name.

   ln -sf $FILE pop.r.nc
   
   echo "pop.r.nc"       >! rpointer.ocn.restart
   echo "RESTART_FMT=nc" >> rpointer.ocn.restart

   echo "Converting restart $FILE ..."

   #----------------------------------------------------------------------
   # Convert each restart file to a (single) DART restart file that has the 
   # timetag for the state ... which we want to overwrite with a generic one.
   # each restart file is called assim_model_state_ud
   #----------------------------------------------------------------------

   ${DARTEXE}/pop_to_dart || exit 1

   #----------------------------------------------------------------------
   # Use the input.nml:restart_file_tool_nml value for the 'valid time'
   # The expected input filename is : dart_input
   # The         output filename is : dart_output
   #----------------------------------------------------------------------
   mv dart_ics dart_input

   ${DARTEXE}/restart_file_tool || exit 2

   #----------------------------------------------------------------------
   # Rename 'generic' output file to reflect the ensemble member number.
   # At the end of the day, we need three things:
   # 1) a DART ics file for each ensemble member :  ens_mem.xxxx 
   # 2) a pointer file  for each ensmeble member :  rpointer.ocn.xxx.restart
   # 3) a POP netCDF restart file (for each ensemble member) actually,
   #    since we're just grabbing the grid dimensions and calendar info
   #    we could probably use just one, but we've got them, and we're just 
   #    creating links ... it's easy enough.
   #----------------------------------------------------------------------
   @ member = $member + 1
   set OFNAME = `printf ens_mem.%03d $member`

   mv dart_output $OFNAME
   echo $FILE:t          >! rpointer.ocn.${member}.restart
   echo "RESTART_FMT=nc" >> rpointer.ocn.${member}.restart
   ln -sf $FILE $FILE:t

   #----------------------------------------------------------------------
   # remove input files to prep for next iteration
   #----------------------------------------------------------------------

   \rm -f pop.r.nc dart_input

end

#-------------------------------------------------------------------------
# Convert restart files to a (single) DART restart file.
#-------------------------------------------------------------------------

# cat ens_mem_* >! filter_ics

#-------------------------------------------------------------------------
# Clean out the big files that are no longer needed.
#-------------------------------------------------------------------------

#\rm -f ens_mem_*

exit 0


