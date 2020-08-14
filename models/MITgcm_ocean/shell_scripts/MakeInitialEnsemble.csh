#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# So the background is that I want 40 snapshot files to create an initial
# conditions file with 40 ensemble members. 
#
# For lack of a better idea: I ran the model starting from 1996 01 01 ??Z 
# for 20 days with a 900 second timestep and a dumpFreq of 12 hours.
# This actually generates 41 snapshot files, but who's counting.
# I just deleted the LAST snapshot before running this script.
#
# We're just going to make the first snapshot into ensemble member 1. 
# Second snapshot will be ensemble member 2 ... I'm just not going to use
# the last snapshots. This may not be the greatest idea - I don't know
# how much spread we will have to start. There must be SOME ...
# depends on the nonlinearity of the model.
#
# Repeatedly call trans_pv_sv to turn each snapshot into a 
# DART initial condition file. Each DART initial conditions file 
# starts with a timestamp that reflects the valid time for the ensuing 
# model state. To generate an ensemble of initial conditions files
# from successive snapshot files ... we will need to overwrite this
# timestamp so that ALL the DART initial conditions files reflect 
# the same timestamp rather than the true timestamp.
#
# The call to restart_file_utility changes the time tag in the file
# to be whatever is in the namelist. In this case, we want the time
# tag to be 1996 01 01 00Z ... in DART-speak ... 144270 days 0 seconds.
#
# TJH - 28 July 2008

# WARNING: You should run this in an 'empty' directory
# WARNING: You should run this in an 'empty' directory
# WARNING: You should run this in an 'empty' directory

# Get input for this execution

set SNAPSHOTDIR = /ptmp/thoar/MITgcm/nodart/fortnight
set     DARTEXE = /fs/image/home/thoar/SVN/DART/models/MITgcm_ocean/work
set     DARTNML = /fs/image/home/thoar/SVN/DART/models/MITgcm_ocean/work/input.nml

# Check to see if all the pieces exist

set checkstat = 0

if ( -d ${SNAPSHOTDIR} ) then
   echo "Using snapshot files from ${SNAPSHOTDIR}"
else
   echo "ERROR:snapshot files directory ${SNAPSHOTDIR} does not exist."
   set checkstat = 1
endif

if ( -f ${SNAPSHOTDIR}/data ) then
   echo "Using ${SNAPSHOTDIR}/data"
else
   echo "ERROR: ${SNAPSHOTDIR}/data does not exist."
   set checkstat = 1
endif

if ( -f ${SNAPSHOTDIR}/data.cal ) then
   echo "Using ${SNAPSHOTDIR}/data.cal"
else
   echo "ERROR: ${SNAPSHOTDIR}/data.cal does not exist."
   set checkstat = 1
endif

if ( -f ${DARTNML} ) then
   echo "Using ${DARTNML}"
else
   echo "ERROR: ${DARTNML} does not exist."
   set checkstat = 1
endif

if ( -f ${DARTEXE}/trans_pv_sv ) then
   echo "Using DART trans_pv_sv from ${DARTEXE}"
else
   echo "ERROR: ${DARTEXE}/trans_pv_sv does not exist."
   echo " Build ${DARTEXE}/trans_pv_sv"
   set checkstat = 1
endif

if ( -f ${DARTEXE}/restart_file_utility ) then
   echo "Using DART restart_file_utility from ${DARTEXE}"
else
   echo "ERROR: ${DARTEXE}/restart_file_utility does not exist."
   echo " Build ${DARTEXE}/restart_file_utility"
   set checkstat = 1
endif

if ( $checkstat > 0 ) then
   echo "ERROR: required input failure."
   exit 9
endif

#-------------------------------------------------------------------------
# The whole operation relies on the fact that there is a directory
# somewhere that has the snapshot files and the corresponding input
# files: data, data.cal - as well as (potentially) another directory
# that has the DART namelists and executables.
#
# We will create a little subdirectory off the DARTEXE directory to
# contain the DART restart files, etc.
#-------------------------------------------------------------------------
# trans_pv_sv needs three namelist files:  data, data.cal and input.nml

cp -p ${SNAPSHOTDIR}/data        .
cp -p ${SNAPSHOTDIR}/data.cal    .
cp -p ${DARTNML}                 .

#-------------------------------------------------------------------------
# Ensure the namelists have the right values.
#-------------------------------------------------------------------------
# We need to run the editor in batch mode.  If you have 'vim' it needs
# one flag; if you have only the older vanilla 'vi' you need another.
# On several systems 'vi' is a link to 'vim' and uses the newer syntax
# so you cannot distinguish which flag will be needed based only on name.
# First try to run 'vim' by full name and then back off to plain 'vi'
# if it is not found.  Punt if neither is found.
#-------------------------------------------------------------------------
set VI_EXE = `which vim`
if ( -x "${VI_EXE}" ) then
   setenv VI 'vim -e'
else
   set VI_EXE = `which vi`
   if ( -x "${VI_EXE}" ) then
      setenv VI 'vi -s'
   else
      echo ""
      echo "Neither the vim nor the vi editor were found.  This script"
      echo "cannot continue unless it can use one of them to update"
      echo "the test input namelist files."
      echo ""
      exit 2
   endif
endif

# Need to modify rest of input.nml for test run
# Essentially, we want the namelist block to look like:
#&restart_file_utility_nml
#   input_file_name              = "assim_model_state_ud",
#   output_file_name             = "filter_updated_restart",
#   ens_size                     =  1,
#   single_restart_file_in       = .true.,
#   single_restart_file_out      = .true.,
#   write_binary_restart_files   = .true.,
#   overwrite_data_time          = .true.,
#   new_data_days                =  144270,
#   new_data_secs                =       0,
#   input_is_model_advance_file  = .false.,
#   output_is_model_advance_file = .false.,
#   overwrite_advance_time       = .false.,
#   new_advance_days             =  -1,
#   new_advance_secs             =  -1   /

echo ':0'                                      >! vi_script
echo '/restart_file_utility_nml'               >> vi_script
echo '/input_file_name'                        >> vi_script
echo ':s/filter_restart/assim_model_state_ud/' >> vi_script
echo '/write_binary_restart_files'             >> vi_script
echo ':s/.false./.true./'                      >> vi_script
echo '/overwrite_data_time'                    >> vi_script
echo ':s/.false./.true./'                      >> vi_script
echo '/new_data_days'                          >> vi_script
echo ':s/-1/144270/'                           >> vi_script
echo '/new_data_secs'                          >> vi_script
echo ':s/-1/0/'                                >> vi_script
echo ':wq'                                     >> vi_script

( ${VI} input.nml < vi_script)

\rm -f vi_script

#-------------------------------------------------------------------------
# Loop over all the snapshot files.
#-------------------------------------------------------------------------

@ memcount = 0

foreach FILE ( ${SNAPSHOTDIR}/T.*.data )

   set FILENAME = $FILE:t
   set FILEBASE = $FILENAME:r
   set TIMESTEP = $FILEBASE:e

   ln -s ${SNAPSHOTDIR}/*.${TIMESTEP}.*   .

   echo "Converting snapshot timestep $TIMESTEP ..."

   #----------------------------------------------------------------------
   # Convert each snapshot file to a (single) DART restart file that has the 
   # timetag for the state ... which we want to overwrite with a generic one.
   #----------------------------------------------------------------------

   echo $TIMESTEP | ${DARTEXE}/trans_pv_sv

   #----------------------------------------------------------------------
   # Use the input.nml:restart_file_utility_nml value for the 'valid time'
   #----------------------------------------------------------------------

   ${DARTEXE}/restart_file_utility 

   #----------------------------------------------------------------------
   # Rename 'generic' output file to reflect the ensemble member number.
   #----------------------------------------------------------------------
   @ memcount = $memcount + 1

   if ( $memcount < 10 ) then
      set OFNAME = ens_mem_00$memcount
   else if ( $memcount < 100 ) then
      set OFNAME = ens_mem_0$memcount
   else
      set OFNAME = ens_mem_$memcount
   endif

   mv filter_updated_restart $OFNAME

   #----------------------------------------------------------------------
   # remove input files to prep for next iteration
   #----------------------------------------------------------------------

   \rm -f *.${TIMESTEP}.* assim_model_state_ud

end

#-------------------------------------------------------------------------
# Convert snapshot files to a (single) DART restart file.
#-------------------------------------------------------------------------

cat ens_mem_* >! filter_ics

#-------------------------------------------------------------------------
# Clean out the big files that are no longer needed.
#-------------------------------------------------------------------------

\rm -f ens_mem_*

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

