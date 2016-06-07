#!/usr/bin/env python
#
# DART software - Copyright 2004 - 2015 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# CREDIT: This script was donated to DART by Luke Madaus during his time
# at the University of Washington. Thanks Luke!
#
# DART $Id$
#
# Modified for Python (Oct. 2015) Luke Madaus, University of Washington
#-----------------------------------------------------------------------------
# Script for use in assimilation applications
# where the model advance is executed as a separate process.
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
# 1) determine how many ensemble members (num_states) this process 
#    will need to advance.
# 2) convey the forecast length to the model
# 3) run the model (make the forecast)
# 4) determine if there are more ensemble members to advance

from __future__ import print_function, division
import os, sys
from netCDF4 import Dataset
from datetime import datetime, timedelta
from namelist_utils import read_namelist, write_namelist

# Read from the input arguments
myname, process, num_states, control_file = sys.argv
process = int(process)
num_states = int(num_states)

print("advance_model: process {:d}".format(process))
print("advance_model: number of states to advance {:d}".format(num_states))
print("advance_model: control file name is {:s}".format(control_file))

centraldir = os.getcwd() 

# Loop through each model state
# Python indexing starts at 0
ensemble_member_line = 0 
input_file_line = 1
output_file_line = 2
current_time_line = 3
advance_to_time_line = 4

state_copy = 1

while state_copy <= num_states:

    #-------------------------------------------------------------------
    # Block 1: Parse the information from the control file
    #
    # NOTE: the input_file and output_file are not used by this script
    #       because of the implementation of the input_filelist.txt constuct.
    #       They are here merely to maintain consistency with old versions
    #       of advance_model.csh scripts. Could easily be removed.
    #-------------------------------------------------------------------

    #set ensemble_member = `head -n $ensemble_member_line  $control_file | tail -n 1`
    #set input_file      = `head -n $input_file_line       $control_file | tail -n 1`
    #set output_file     = `head -n $output_file_line      $control_file | tail -n 1`
    #set current_time    = `head -n $current_time_line     $control_file | tail -n 1`
    #set advance_to_time = `head -n $advance_to_time_line  $control_file | tail -n 1`
    with open(control_file,'r') as cntrl:
        cntrl_lines = cntrl.readlines()
        # Strip the \n off the ends...
        ensemble_member = int(cntrl_lines[ensemble_member_line][:-1])
        input_file      = cntrl_lines[input_file_line][:-1]
        output_file     = cntrl_lines[output_file_line][:-1]
        current_time    = cntrl_lines[current_time_line][:-1]
        advance_to_time = cntrl_lines[advance_to_time_line][:-1]





    # The run-time directory for the entire experiment is called CENTRALDIR;
    # we need to provide a unique directory for each model advance.
    temp_dir = 'dir_model_{:03d}'.format(ensemble_member)
    # Don't know if we need this yet...probably don't...
    if not os.path.exists(temp_dir):
        os.system('mkdir {:s}'.format(temp_dir))
    os.chdir(temp_dir)
    os.system('cp ../namelist.input namelist.input.template'.format(temp_dir))
    os.system('ln -sf ../LANDUSE.TBL .')
    os.system('ln -sf ../input_sounding .')

    #-------------------------------------------------------------------
    # Block 2:  Convey the forecast length to the model.
    #
    # Create two time strings that we can difference to get the number
    # of seconds to advance which must be conveyed to CM1. We need to
    # parse the times from the control_file (whose format is determined
    # by DART in consideration of all the models) to a time useful for CM1
    #-------------------------------------------------------------------

    # This whole thing is redone with python datetime
    currentdt = datetime(*map(int, current_time.split()[1:]))
    advancedt = datetime(*map(int, advance_to_time.split()[1:]))
    # Compute the difference
    diffdt = advancedt - currentdt
    totaldays = diffdt.days
    totalsecs = diffdt.seconds
    run_time = totaldays * 86400 + totalsecs


    """
    # The dart function 'advance_time' needs an input.nml
    os.system('ln -svf  ../input.nml  .  || exit 1')
    set bob=`echo $current_time`
    set year   = $bob[2]
    set month  = $bob[3]
    set day    = $bob[4]
    set hour   = $bob[5]
    set minute = $bob[6]
    set second = $bob[7]
    set cur_datestr = `printf %04d%02d%02d%02d%02d%02d $year $month $day $hour $minute $second`
    set current_d=(`echo ${cur_datestr} 0 -g | ../advance_time`)

    #echo "We think the current  date string is: ${cur_datestr}"
    #echo "Current  days and seconds: $current_d"
    #echo "Started with: ${current_time}"

    set bob=`echo $advance_to_time`
    set year   = $bob[2]
    set month  = $bob[3]
    set day    = $bob[4]
    set hour   = $bob[5]
    set minute = $bob[6]
    set second = $bob[7]
    set fcst_datestr = `printf %04d%02d%02d%02d%02d%02d $year $month $day $hour $minute $second`
    set forecast_d =(`echo ${fcst_datestr} 0 -g | ../advance_time`)

    #echo "We think the forecast date string is: ${fcst_datestr}"
    #echo "Forecast days and seconds: $forecast_d"
    #echo "Started with: ${advance_to_time}"

    @ totaldays = ( $forecast_d[1] - $current_d[1] )
    @ totalsecs = ( $forecast_d[2] - $current_d[2] )

    @ run_time = $totaldays * 86400 + $totalsecs
    """


    #echo "DEBUG:Run time: $run_time"

    # Now that we know the forecast length, must convey that to CM1 by
    # replacing the bogus string with the real forecast length.
    # Both run_time and rstfrq (restart frequency) are changed.
    # FIXME: check the impact of changing the restart frequency
    # has on whole thing by varying the forecast lengths.
    #
    # NOTE: irst == 1, always. Indicates using an existing state.
    # By forcing rstnum == 1, we can always link the most current
    # file to be cm1out_rst_000001.nc
    
    # Another adjustment with python, now using namelist_utils
    os.system('rm -f namelist.input')
    nmld = read_namelist('namelist.input.template')
    # Change the foreacst and restart length
    nmld['param1']['run_time'] = run_time
    nmld['param1']['rstfrq'] = run_time
    # Also, for serial run, make sure nodex and nodex are 1
    nmld['param0']['nodex'] = 1
    nmld['param0']['nodey'] = 1
    # Write the new namelist
    write_namelist(nmld, 'namelist.input')


    """
    rm -f namelist.input
    sed -e "s/CM1_FORECAST_LENGTH/${run_time}/" namelist.input.template > namelist.input

    grep CM1_FORECAST_LENGTH namelist.input
    if ($status == 0) then
      echo "The CM1 namelist file 'namelist.input' did not get the new run_time."
      echo "Aborting."
      exit 2
    endif
    """

    # LEM Here is where I'm putting the initial file turnaround
    # For the first time, look for the restart file in the control
    # directory that matches this ensemble member at currentdt
    possible_first_file = 'cm1out.{:03d}.{:%Y%m%d%H%M%S}.nc'.format(ensemble_member, currentdt)
    if not os.path.exists('cm1out_rst_000001.nc') and os.path.exists('../{:s}'.format(possible_first_file)):
        os.system('cp ../{:s} cm1out_rst_000001.nc'.format(possible_first_file))
    


    #-------------------------------------------------------------------
    # Block 3: Advance the model.
    #
    # Saving the run-time messages from each file to a unique log file.
    # This is intended to make debugging easier.
    #
    # CM1 always expects a restart file called "cm1out_rst_000001.nc".
    #
    # DART is going to save the most current state with a date/time tag 
    #    and must also link that most current state to the static 
    #    "cm1out_rst_000001.nc" name.
    #-------------------------------------------------------------------
    fcst_datestr = '{:%Y%m%d%H%M%S}'.format(advancedt)
    os.system('rm -rf cm1log.{:s}.txt'.format(fcst_datestr))

    os.system('../cm1.exe |& tee cm1log.{:s}.txt'.format(fcst_datestr))

    # check the file
    with open('cm1log.{:s}.txt'.format(fcst_datestr),'r') as logfile:
        contents = logfile.read()
        if 'Program terminated normally' not in contents:
            print("ERROR: cm1 model advance failed.")
            print("ERROR: check {:s}/cm1log.{:s}.txt".format(temp_dir, fcst_datestr))
            exit(3)
        else:
            print("CM1 should now be at: {:s}".format(fcst_datestr))

            
    """
    grep 'Program terminated normally' cm1log.${fcst_datestr}.txt
    if ($status != 0) then
      echo "ERROR: cm1 model advance failed."
      echo "ERROR: check $temp_dir/cm1log.${fcst_datestr}.txt"
      exit 3
    else
      echo "CM1 should now be at $fcst_datestr"
    endif
    """

    # Trying to grab the most recent restart file
    # Python os listing here to sort files by order
    outfiles = [f for f in os.listdir('.') if f.startswith('cm1out_rst_') and f.endswith('.nc')]
    outfiles.sort()
    outfile = outfiles[-1]

    # Another python check...just to be sure, open this output file
    # and compute its date to compare with advancedt
    # Read the simulation start (base time) from the namelist
    nmldtime = nmld['param11']
    startdt = datetime(nmldtime['year'], nmldtime['month'], nmldtime['day'],\
                   nmldtime['hour'], nmldtime['minute'], nmldtime['second'])
    # Now get the time from what we think is the latest restart file
    with Dataset(outfile, 'r') as rst:
        filesec = int(rst.variables['time'][0])
    # See if this matches where we think we should be
    filedt = startdt + timedelta(seconds=filesec)
    if filedt != advancedt:
        print("ERROR: cm1 model advance not at right time")
        print("ERROR: we think this is latest file: {:s}/{:s}".format(temp_dir, outfile))
        print("ERROR: but filetime is {:%Y%m%d%H%M%S} and desired time is {:%Y%m%d%H%M%S}".format(filedt, advancedt))
        exit(3)


    # Rename this to the "official" file
    os.system('mv -v {:s} cm1out.{:03d}.{:s}.nc'.format(outfile, ensemble_member, fcst_datestr))

    # Update filename referenced by CENTRALDIR/input_filelist.txt
    # to point to the most recent restart file.
    os.system('ln -svf cm1out.{:03d}.{:s}.nc cm1out_rst_000001.nc'.format(ensemble_member, fcst_datestr))

    #-------------------------------------------------------------------
    # Block 4:
    # Update the location in the control file with the information
    # for the next ensemble member
    #-------------------------------------------------------------------

    state_copy           += 1
    ensemble_member_line += 5
    input_file_line      += 5
    output_file_line     += 5
    current_time_line    += 5
    advance_to_time_line += 5

    # Change back to original directory
    os.chdir(centraldir)

# Return to normal indent---end of while block

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

os.system('rm -rf {:s}'.format(control_file))

exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

