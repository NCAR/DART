#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.
# Can be used as-is with most low-order models and the bgrid model which
# can be advanced using the integrate_model executable.
# 
# Arguments are (created by 'filter' or 'perfect_model_obs' and include):
# 1) the process number of caller,
# 2) the number of ensemble members/state copies belonging to that process, and 
# 3) the name of the control_file for that process.
# 
# If this script finishes and the 'control_file' still exists, it is
# an ERROR CONDITION and means one or more of the ensemble members did
# not advance properly. Despite our best attempts to trap on this
# condition, some MPI installations simply hang, some properly terminate.
#
# This script loops over all the entries in the control_file to advance 
# any/all of the ensemble members.  The number of trips through the 
# loop is the second argument to this script. The control_file contains 
# the information about which ensemble members are to be advanced by THIS TASK.
# Sometimes it may be just one ensemble member, sometimes all of them.
# Read DART/doc/html/filter_async_modes.html and the mpi_intro.html
# for an overview.
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance
#    and copies the necessary files into the temporary directory
# 2) copies/converts the DART state vector to something the model can ingest
# 3) runs the model
# 4) copies/converts the model output to input expected by DART

set      process = $1
set   num_states = $2
set control_file = $3

#----------------------------------------------------------------------
# Block 1: copy necessary input files/executables/files common
#          to all model advances to a clean, temporary directory.
#          These will be used by ALL of the ensemble
#          members being advanced by this script.
#----------------------------------------------------------------------

# Create a unique temporary working directory for this process's stuff
# The run-time directory for the entire experiment is called CENTRALDIR;
# we need to provide a safe haven for each TASK ... in 'temp_dir'.

set temp_dir = 'advance_temp'$process

# Create a clean temporary directory and go there
# \rm -rf  $temp_dir  || exit 1
mkdir -p $temp_dir  || exit 1
cd       $temp_dir  || exit 1

# Get the DART input.nml and the NOAH namelist

foreach FILE ( GENPARM.TBL SOILPARM.TBL URBPARM.TBL VEGPARM.TBL namelist.hrldas input.nml )
   \cp -v ../$FILE . || exit 2
end

set MYSTRING = `grep HRLDAS_CONSTANTS_FILE namelist.hrldas`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set NOAHFILE = `echo $MYSTRING[$#MYSTRING]`
\ln -sv ../$NOAHFILE .

# Extract the directory containing the forcing files used to create the LDASIN files.
# This script uses logic that requires ALL of the hourly forcing files
# to be concatenated into a single netCDF file that uses 'time' as the unlimited
# dimension. The required time slices are extracted from this single netCDF file
# into the filenames expected by NOAH. Since each ensemble member generally gets
# a unique atmospheric forcing, this helps minimize the number of files. 

set MYSTRING   = `grep FORCING_FILE_DIRECTORY namelist.hrldas`
set MYSTRING   = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
set FORCINGDIR = `echo $MYSTRING[$#MYSTRING]`

# echo "Master atmospheric netCDF forcing file(s) coming from $FORCINGDIR"

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`
   set instance        = `printf "%04d" $ensemble_member`

   #-------------------------------------------------------------------
   # Block 2: copy/convert the DART state vector to something the 
   #          model can ingest.
   #
   #          * remove any NOAH scraps from previous advances
   #          * copy/link ensemble-member-specific files
   #          * convey the advance-to-time to the model
   #          * convert the DART state vector to model format 
   #-------------------------------------------------------------------

   # echo "advance_model.csh block 2 converting ensemble member $instance"

   \rm -f *.LDASIN* *.LDASOUT*
   \rm -f restart.nc dart_restart noah_advance_information.txt
   \ln -v ../restart.$instance.nc  restart.nc   || exit 2
   \ln -s ../$input_file           dart_restart || exit 2
   ../dart_to_noah                              || exit 2

   if ( ! -e noah_advance_information.txt ) then
      echo "ERROR: dart_to_noah failed for member $ensemble_member"
      echo "ERROR: dart_to_noah failed for member $ensemble_member"
      echo "ERROR: dart_to_noah failed for member $ensemble_member"
   endif

   # This next two parts are based on using one-hour forcing files
   # since the minimum time to advance the model seems to be 1 hour.
   # (kday, khour, but no kminute, for example)
   # dart_to_noah provides the setting for namelist.hrldas:khour
   # we need to put that value in the local copy of namelist.hrldas

   set numadvancestr = `grep -i khour noah_advance_information.txt`
   set numadvancestr = `echo $numadvancestr | sed -e "s#[=,']# #g"`
   set numadvancestr = `echo $numadvancestr | sed -e 's#"# #g'`
   set numadvances   = `echo $numadvancestr[$#numadvancestr]`

ex namelist.hrldas <<ex_end
g;KHOUR ;s;= .*;= $numadvances;
wq
ex_end

   # The forcing has to be for the NEXT "FORCING_TIMESTEP", apparently.
   # FORCING_TIMESTEP is defined in namelist.input At this point, dart_to_noah
   # has assumptions that the forcing_timestep is one hour.

   set numfilestring = `grep -ni nfiles noah_advance_information.txt`
   set numfilestring = `echo $numfilestring | sed -e "s#[=,':]# #g"`
   set numfilestring = `echo $numfilestring | sed -e 's#"# #g'`
   set numfiles      = `echo $numfilestring[$#numfilestring]`
   set skipNlines    = `echo $numfilestring[1]`

   # Extract needed forcing periods from master NetCDF file located in $FORCINGDIR.
   #
   # For PERTURBED FORCING ... there must be a forcing file for each ensemble member
   # with names like 'noah_forcing.$FYEAR.nnnn.nc' - where nnnn is a 
   # zero-filled integer corresponding to the ensemble member instance [0001, 0002, ...]

   @ ifile = 1
   while ($ifile <= $numfiles)
      @ linenum = $skipNlines + $ifile
      set FNAME = `head -$linenum noah_advance_information.txt | tail -1`
      set FDATE = `echo $FNAME | sed -e "s#[.,']# #g"`
      set FDATE = `echo $FDATE[1]`
      set FYEAR = `echo $FDATE] | cut -c1-4`
      set FFILE = $FORCINGDIR/noah_forcing_$FYEAR.$instance.nc
   
      # Print some message
      # echo 'extracting forcing for ' $FDATE.LDASIN_DOMAIN1

      # Now create forcing data for single timestep
      ncks -O -a -d time,$FDATE. $FFILE $FDATE.LDASIN_DOMAIN1

      if ($status != 0) then
         echo "ERROR: cannot create LDASIN file for $FDATE from $FFILE"
         echo "ERROR: cannot create LDASIN file for $FDATE from $FFILE"
         exit 20
      endif

      @ ifile = $ifile + 1
   end

   #-------------------------------------------------------------------
   # Block 3: advance the model
   #          In this case, we are saving the run-time messages to
   #          a LOCAL file, which makes debugging easier.
   #          integrate_model is hardcoded to expect input in temp_ic 
   #          and it creates temp_ud as output. 
   #          Your model will likely be different.
   #-------------------------------------------------------------------

   ../Noah_hrldas_beta

   set noah_status = `ls -1 RESTART*DOMAIN* | wc -l`
   if ($noah_status < 1)  then
      echo "ERROR: NOAH died"
      echo "ERROR: NOAH died"
      ls -l
      exit 23
   else if ($noah_status > 1) then
      echo "WARNING: NOAH created the following RESTART files. only expected one." 
      ls -l RESTART*DOMAIN*
   endif

   \rm -f restart.nc

   #-------------------------------------------------------------------
   # Block 4: convert the new model state into a DART-readable form. 
   #          rename the restart file to reflect the ensemble member ID
   #
   # We want the LAST restart file and the FIRST ldasout file.
   #-------------------------------------------------------------------

   set RESTART = `ls -1  RESTART* | tail -1`
   set LDASOUT = `ls -1 *LDASOUT* | head -1`

   \ln -s $RESTART  restart.nc  || exit 4
   ../noah_to_dart              || exit 4

   \mv -v  dart_ics  ../$output_file          || exit 5
   \mv -v  $RESTART  ../$RESTART.$instance.nc || exit 5
   \mv -v  $LDASOUT  ../$LDASOUT.$instance.nc || exit 5
   \ln -sf $RESTART.$instance.nc ../restart.$instance.nc
   \ln -sf $RESTART.$instance.nc ../restart.nc

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3

end

# Change back to original directory and get rid of temporary directory.
# If all goes well, there should be no need to keep this directory.
# If you are debugging, you may want to keep this directory. 

cd ..
\rm -rf $temp_dir

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

\rm -rf $control_file

exit 0


