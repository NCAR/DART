#!/usr/bin/env python
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# CREDIT: This script was donated to DART by Luke Madaus during his time
# at the University of Washington. Thanks Luke!
#
# Modified for Python (Oct. 2015) Luke Madaus, University of Washington
#-----------------------------------------------------------------------------
#
# This script is meant to provide an example of what is necessary to set up
# a very simple experiment that involves DART and a very small ensemble of
# CM1 model states. This is guaranteed to be scientifically bogus - it is just
# a test of the file motion and mechanics.
#
# ASSUMPTIONS/REQUIREMENTS
#
# *) CM1 is required to use netCDF restart files.
#
# *) A collection of CM1 model states for initial conditions will be available
#
# *) There is a separate observation sequence file for each assimilation time
#
# *) The DART input.nml file has some required values as defined below.
#
# *) Each time CM1 is advanced, it will start from the same filename,
#    and the restart number in that filename will be 000001 - ALWAYS.
#    That filename will be a link to the most current model state.
#
#-----------------------------------------------------------------------------
# These are the DART namelist variables that MUST have these values.
#
# &filter_nml
#   input_state_file_list       = "input_filelist.txt"
#   output_state_file_list      = "output_filelist.txt"
#   init_time_days           = -1
#   init_time_seconds        = -1
#   obs_seq_in_file_name     = "obs_seq.out"
#   ...
# /
# 
#-----------------------------------------------------------------------------
# Python imports here
import os
import argparse
from netCDF4 import Dataset
from datetime import datetime, timedelta
from namelist_utils import write_namelist, read_namelist

# There are two command line arguments
parser = argparse.ArgumentParser(description='Setup filter for given model time')

# This tells filter what the current assimilation time should be
parser.add_argument('-d', dest='filtertime', type=str, help='Date of assimilation in YYYYMMDDHHMMSS')

# This flag tells the script whether or not we are
# initializing the ensemble experiment.  Default (not using this flag)
# is to assume the ensemble is already set up

parser.add_argument('-i', dest='initialize', action='store_true', default=False,\
    help="Flag to indicate that the script should search for ensemble members in external dir.")

args = parser.parse_args()
try:
    assim_target_date = datetime.strptime(args.filtertime, '%Y%m%d%H%M%S')
except:
    print("An analysis time is required to setup filter.")
    print("Please run setup_filter.py with the flag:")
    print("-d YYYYMMDDHHmmss for the current analysis time")
    exit(1)
# Are we initializing from remote directory?
initialize = args.initialize

###############################################
#       BEGIN USER-DEFINED VARIABLES          #
###############################################
# Set some parameters about our experiment
jobname            = 'cm1_filter'             # Name for this experiment
ens_size           = 3                        # Number of ensemble members
restart_filename   = 'cm1out_rst_000001.nc'   # Restart file name
window_mins        = 15                       # Assimilation window (+/-minutes)

# some systems don't like the -v option to any of the following

copy    = 'cp -pv'
link    = 'ln -sv'
remove  = 'rm -rf'

# Information about the queue we are submitting to
queue_system = 'BSUB'
mpi_run_command = 'mpirun.lsf'
queue_sub_command = 'bsub'
jobsub_info = {'-P' : 'P8685005x',
              '-W' : '00:10',
              '-n' : '16',
              '-q' : 'small',
              '-J' : 'run_cm1',
              '-o' : 'run_cm1.%J.out',
              '-e' : 'run_cm1.%J.err'}

# List of DART output files to archive after each cycle
files_to_archive = ['obs_seq.final', \
   'preassim_mean.nc', 'preassim_sd.nc','preassim_priorinf_mean.nc', 'preassim_priorinf_sd.nc', \
   'analysis_mean.nc', 'analysis_sd.nc','analysis_priorinf_mean.nc', 'analysis_priorinf_sd.nc', \
     'output_mean.nc',   'output_sd.nc',  'output_priorinf_mean.nc',   'output_priorinf_sd.nc']

#-----------------------------------------------------------------------------
# Define the directory that will hold the experiment (i.e. centraldir).
# All the necessary resources get staged here.
# All the run-time input should be defined/modified here.
# This directory will be overwritten if this script is run in initialize mode
#-----------------------------------------------------------------------------

#centraldir = '{home}/temp/{JOBNAME}/job_{JOBID}'.format(**locals())
centraldir = '/glade/scratch/thoar/cm1_test/{:s}'.format(jobname)

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
# dartdir      The location of the DART cm1 model directory
# cm1dir       The location of the cm1 executable
# icdir        The location of the initial ensemble of cm1 files (if initialize==True)
# obsdir       The location of the observation sequence files
#-----------------------------------------------------------------------------

dartdir = '/glade/p/work/thoar/DART/rma_trunk/models/cm1'
cm1dir = '/glade/u/home/lmadaus/cm1/cm1r18/run'
icdir = '/glade/p/work/lmadaus/cm1/cm1_example/mems'
obsdir = '/glade/p/work/lmadaus/cm1/cm1_example/obs'

###############################################
#         END USER-DEFINED VARIABLES          #
###############################################

# This is an odd pythonic thing for joining a list of strings into a single string
archive_filestr = "','".join(files_to_archive)

# Make the working directory
print(jobname, 'centraldir ==', centraldir)

if initialize:
    print("INITIALIZING NEW WORK DIRECTORY")
    if os.path.exists(centraldir):
        os.system('{remove} {centraldir}'.format(**locals()))
    os.system('mkdir -p {:s}'.format(centraldir))
    os.system('mkdir {:s}/archive'.format(centraldir))
os.chdir(centraldir)

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the cm1 executable, control files, and data files.
# Only if initializing
#-----------------------------------------------------------------------------
if initialize:
    os.system('{copy} {dartdir}/work/filter          . || exit 1'.format(**locals()))
    os.system('{copy} {dartdir}/work/advance_time    . || exit 1'.format(**locals()))
    os.system('{copy} {dartdir}/work/input.nml       . || exit 1'.format(**locals()))
    os.system('{copy} {cm1dir}/cm1.exe               . || exit 1'.format(**locals()))

# Find the observation sequence file at the right time
if not os.path.exists('{:s}/{:%Y%m%d%H%M%S}_obs_seq.prior'.format(obsdir, assim_target_date)):
    print("Could not find obs sequence file:")
    print("{:s}/{:%Y%m%d%H%M%S}_obs_seq.prior".format(obsdir, assim_target_date))
    print("Aborting.")
    exit(1)
os.system('{copy} {obsdir}/{assim_target_date:%Y%m%d%H%M%S}_obs_seq.prior obs_seq.prior'.format(**locals()))

# Make a subdirectory in the archive directory to hold output for this cycle
os.system('mkdir {:s}/archive/{:%Y%m%d%H%M%S}'.format(centraldir, assim_target_date))

# Load the input.nml file and change it to reflect this time
dartnml = read_namelist('input.nml')

#> if not initializing, we need to ensure that the output inflation files from
#> previous cycle are used as input for next cycle
if not initialize:
	dartnml['filter_nml']['inf_initial_from_restart'] = [True, True]
	dartnml['filter_nml']['inf_sd_initial_from_restart'] = [True, True]

# Change the name of the template file
dartnml['model_nml']['cm1_template_file'] = restart_filename

# Set the right obs sequence name
dartnml['filter_nml']['obs_sequence_in_name'] = 'obs_seq.prior'

# Set the number of ensemble members
dartnml['filter_nml']['ens_size'] = ens_size

# Set some interesting observation types to assimilate
dartnml['obs_kind_nml']['assimilate_these_obs_types'] = ['TEMPERATURE_2M']
#>@todo TJH test the rest of the observation types in the obs file
#dartnml['obs_kind_nml']['evaluate_these_obs_types'] = ['TEMPERATURE_2M']

# Write this modified namelist out
write_namelist(dartnml, 'input.nml')

# filter will read the contents of input_filelist.txt for the names
# of the CM1 restart files. After each model advance these same
# filenames will be queried for the new model state.
# So, the actual file associated with this name must be updated.
# Ensemble advances are handled by running advance_ensemble.py
# Each model advance is performed in a unique model directory that
# has an instance or ensemble member number in the directory name.

os.system('{remove} input_filelist.txt'.format(**locals()))
os.system('touch input_filelist.txt')

for instance in range(1,ens_size+1):

    # filter will read the contents of input_filelist.txt for
    # the file containing the most up-to-date model state.
    # Each ensemble member (instance) will get advanced in its own directory.

    dirname = 'dir_model{:03d}'.format(instance)
    # If we are in initialize mode, copy the initial conditions from icdir
    # Otherwise, just make sure that a restart file is in the directory for
    # this member
    if initialize:
        os.system('mkdir {:s}'.format(dirname))
        os.system('{copy} {icdir}/m{instance}/cm1out_rst_000001.nc {dirname}/{restart_filename} || exit 2'.format(**locals()))
        os.system('{copy} {icdir}/m{instance}/input_sounding       {dirname} || exit 2'.format(**locals()))
        os.system('{copy} {icdir}/m{instance}/LANDUSE.TBL          {dirname} || exit 2'.format(**locals()))
        os.system('{copy} {icdir}/m{instance}/namelist.input       {dirname}/namelist.input.template || exit 2'.format(**locals()))
        os.system('{link} {cm1dir}/cm1.exe              	   {dirname}/cm1.exe || exit 1'.format(**locals()))
    elif not os.path.exists('{dirname}/{restart_filename}'.format(**locals())):
            print("ERROR: unable to find restart file: {restart_filename} for member {instance}".format(**locals()))
            exit(1)

    
    # If we are here, then there is at least one restart file in the directory
    # Find the restart file at the correct time and name it appropriately
    # Temporarily save the previous restart file as prev_{restart_filename}
    # List all restart files
    rst_files = [f for f in os.listdir(dirname) if f.startswith('cm1out_rst')]
    # Start with the last restart file
    rst_files.sort()
    rst_files.reverse()
    rst_file_found = False
    for rst_file in rst_files:
        if not rst_file_found:
            # Try this file
            with Dataset('{dirname}/{rst_file}'.format(**locals()), 'r') as dset:
                curtime = dset.variables['time'][-1]
                init = datetime(dset.year, dset.month, dset.day, dset.hour, dset.minute, dset.second)
                validtime = init + timedelta(seconds=int(curtime))
            # Do we match the filter target date?
            if validtime != assim_target_date:
                # Skip to the next file
                continue
            # If we are here, this restart file's time matches the assim_target_date
            # Copy the old restart file (if it exists) to prev_{restart_filename}
            print("Found usable restart file for member {:d}".format(instance))
            if os.path.exists('{dirname}/{restart_filename}'.format(**locals())):
                os.system('{copy} {dirname}/{restart_filename} {dirname}/prev_{restart_filename}'.format(**locals()))
            # Move the newly-found, appropriately-timed restart file to {restart_filename}
            os.system('mv {dirname}/{rst_file} {dirname}/{restart_filename}'.format(**locals()))
            rst_file_found = True
                

    # We've looped through all restart files at this point.  Check to be sure we found one
    # that works
    if not rst_file_found:
        print("No restart files found for member {:d} at requested time!".format(instance))
        print("Requested assimilation date is given as {:%Y%m%d%H%M%S}".format(assim_target_date))
        print("Ensure that the model is at the right time! Aborting")
        exit(1)


    # Append the name of the restart file into the list
    with open('input_filelist.txt','a') as outfile:
        outfile.write('{:s}/{:s}\n'.format(dirname, restart_filename))


# When filter starts, it needs to read the grid information.
# The DART/CM1 interface specifies the file for the grid information by
# the input.nml:&model_nml:cm1_template_file  variable, which defaults
# to cm1out_rst_000001.nc
#
# Consequently, we can get simply link to any model state.
# Until such time as the 'time' variable is correct, we also need a namelist.input
# By having the same files in input_filelist.txt and output_filelist.txt, DART will
# overwrite the file with the posterior from the assimilation. 
os.system('{link} dir_model001/cm1out_rst_000001.nc  cm1out_rst_000001.nc || exit 3'.format(**locals()))
os.system('{copy} dir_model001/namelist.input.template     namelist.input || exit 3'.format(**locals()))
os.system('{copy} input_filelist.txt                  output_filelist.txt || exit 3'.format(**locals()))

# Write a queue submission script for filter
with open('run_filter.py','w') as outfile:
    outfile.write('#!/usr/bin/env python\n')
    # Pythonic way to format all the queue submission info
    for flag, flagval in jobsub_info.items():
        outfile.write("#{:s} {:s} {:s}\n".format(queue_system, flag, flagval))
    # Make sure we are in the right directory
    outfile.write('import os\n')
    outfile.write("os.chdir('{:s}')\n".format(centraldir))
    outfile.write("os.system('{:s} ./filter')\n".format(mpi_run_command))

    #>@todo FIXME check for successful completion some way

    # stage the output inflation files as input for the next assimilation cycle
    outfile.write("\n# staging inflation files for next time if needed \n")
    outfile.write("os.system('{copy} output_priorinf_mean.nc  input_priorinf_mean.nc')\n".format(**locals()))
    outfile.write("os.system('{copy} output_priorinf_sd.nc    input_priorinf_sd.nc  ')\n".format(**locals()))
    outfile.write("os.system('{copy} output_postinf_mean.nc   input_postinf_mean.nc ')\n".format(**locals()))
    outfile.write("os.system('{copy} output_postinf_sd.nc     input_postinf_sd.nc   ')\n\n".format(**locals()))

    # Archive the files if it worked in the queue script
    outfile.write("for f in ['{:s}']:\n".format(archive_filestr))
    outfile.write('    if not os.path.exists(f):\n')
    outfile.write('        print("Error! Could not find DART output file {:s}".format(f))\n')
    outfile.write('        print("Check to be sure DART completed successfully.")\n')
    outfile.write('        continue\n')
    outfile.write('    os.system("mv {:s} archive/{:%Y%m%d%H%M%S}".format(f))\n'.format('{:s}',assim_target_date))

# UNCOMMENT HERE if you want the script to automatically submit this
#os.system('{:s} < run_filter.py'.format(queue_sub_command))

print("")
print("cd {:s}".format(centraldir))
print("Make sure that input.nml contents are correct.")
print("Make sure the obs_seq.in and the model state are appropriate.")
print("Make sure that namelist.input contents are correct.")
print("Make sure that input_filelist.txt contains the right file and the file exists.")
print("Launch (mpirun?) ./filter")
print("Or submit the run_filter.py script in that directory")
print("")

exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

