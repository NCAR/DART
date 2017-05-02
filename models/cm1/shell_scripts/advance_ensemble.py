#!/usr/bin/env python

import os
import argparse
from netCDF4 import Dataset
from datetime import datetime, timedelta
from namelist_utils import write_namelist, read_namelist

# There are two command line argauments
parser = argparse.ArgumentParser(description='Run ensemble members for given integration time')
# This tells filter what the current assimilation time should be
parser.add_argument('-d', dest='initdate', type=str, help='Date of initialization in YYYYMMDDHHMMSS')
parser.add_argument('-l', dest='run_length', type=int, default=3600,\
    help="Length of model integration in seconds")

args = parser.parse_args()
try:
    initdate = datetime.strptime(args.initdate, '%Y%m%d%H%M%S')
except:
    print("An initialization time is required to coordinate runs")
    print("Please run advance_model with the flag:")
    print("-d YYYYMMDDHHmmss for the current initialization time")
    exit(1)
# Get the run length from the argument
run_length = int(args.run_length)

################################################
#       BEGIN USER-DEFINED VARIABLES           #
################################################


# Set some parameters about our experiment--this should match what's in setup_filter.py
jobname            = 'cm1_filter'             # Name for this experiment
ens_size           = 3                        # Number of ensemble members
restart_filename   = 'cm1out_rst_000001.nc'   # Restart file name

# some systems don't like the -v option to any of the following

copy    = 'cp -pv'
link    = 'ln -sv'
remove  = 'rm -rf'

# Information about the queue we are submitting to
queue_system = 'BSUB'
mpi_run_command = 'mpirun.lsf'
queue_sub_command = 'bsub'
jobsub_info = {'-P' : 'P8685005x',
              '-W' : '00:05',
              '-n' : '16',
              '-q' : 'small',
              '-J' : 'run_cm1',
              '-o' : 'run_cm1.%J.out',
              '-e' : 'run_cm1.%J.err'}
              
#-----------------------------------------------------------------------------
# Define the directory that will hold the experiment (i.e. centraldir).
# All the necessary resources get staged here.
# All the run-time input should be defined/modified here.
# This should match what's in setup_filter.py
#-----------------------------------------------------------------------------

#centraldir = '{home}/temp/{JOBNAME}/job_{JOBID}'.format(**locals())
centraldir = '/glade/scratch/thoar/cm1_test/{:s}'.format(jobname)

################################################
#         END USER-DEFINED VARIABLES           #
################################################

print("centraldir = ", centraldir)
# Loop through each ensemble member once to ensure that
# the restart file is there
for member in range(1,ens_size+1):
    print("Prepping member {:d}".format(member))
    memdir = '{:s}/dir_model{:03d}'.format(centraldir, member)
    # Check that the restart file is there and is at the right time
    try:
        with Dataset('{memdir}/{restart_filename}'.format(**locals()), 'r') as dset:
            curtime = dset.variables['time'][-1]
            init = datetime(dset.year, dset.month, dset.day, dset.hour, dset.minute, dset.second)
            validtime = init + timedelta(seconds=int(curtime))
            # Do we match the filter target date?
            if validtime != initdate:
                print("Restart file for member {:d} is at date {:%Y%m%d%H%M%S}".format(member, validtime))
                print("Requested initialization date is given as {:%Y%m%d%H%M%S}".format(initdate))
                print("Ensure that the model is at the right time! Aborting")
                exit(1)
    except:
        # Couldn't open the file
        print("Unable to open file {memdir}/{restart_filename}")
        print("Ensure that there is a model restart file for this member")
        exit(1)

    # Re-write the cm1 namelist for this integration
    cm1nml = read_namelist('{memdir}/namelist.input.template'.format(**locals()))
    # Set the run_length and be sure it knows to start from a restart
    cm1nml['param2']['irst'] = 1
    cm1nml['param2']['rstnum'] = 1
    cm1nml['param1']['run_time'] = run_length
    cm1nml['param1']['rstfrq'] = run_length
    # Write out the updated namelist
    write_namelist(cm1nml, '{memdir}/namelist.input'.format(**locals()))

    # Make a submission file for this member
    with open('{:s}/run_mem.py'.format(memdir), 'w') as outfile:
        outfile.write("#!/usr/bin/env python\n")
        # Pythonic way to format all the queue submission info
        for flag, flagval in jobsub_info.items():
            outfile.write("#{:s} {:s} {:s}\n".format(queue_system, flag, flagval))
        # Make sure we are in the right directory
        outfile.write('import os\n')
        outfile.write("os.chdir('{:s}')\n".format(memdir))
        outfile.write("os.system('{:s} ./cm1.exe')\n".format(mpi_run_command))
        
# Now that we know all members are at the right time, their
# namelists have been updated, and they have a queue submission
# script, submit each member to the queue
this_dir = os.getcwd()
for member in range(1, ens_size+1):
    memdir = '{:s}/dir_model{:03d}'.format(centraldir,member)
    os.chdir(memdir)
    print("Submitting member {:d}".format(member))
    os.system('{:s} < run_mem.py'.format(queue_sub_command))
    os.chdir(this_dir)
 
