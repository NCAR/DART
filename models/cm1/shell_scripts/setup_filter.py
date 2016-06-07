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
#
# This script is meant to provide an example of what is necessary to set up
# a very simple experiment that involves DART and a very small ensemble of
# CM1 model states. This is guaranteed to be scientifically bogus - it is just
# a test of the file motion and mechanics.
#
# ASSUMPTIONS/REQUIREMENTS
#
# *) CM1 is required to be compiled *without* MPI so we can run in async==2
#
# *) CM1 is required to use netCDF restart files.
#
# *) A collection of CM1 model states will be available. Each of those
#    model states is required to have its own namelist.input, LANDUSE.TBL,
#    etc. Each set of states and resources will be in a separate directory
#    named for each ensemble member ... i.e.  xxx/m1/* xxx/m2/* xxx/m3/* ...
#
# *) The DART obs_seq.out file will totally define the duration
#    of the assimilation.
#
# *) The DART input.nml file has some required values as defined below.
#
# *) Each time CM1 is advanced, it will start from the same filename,
#    and the restart number in that filename will be 000001 - ALWAYS.
#    That filename will be a link to the most current model state.
#
# *) Each time CM1 has finished advancing the most current restart file
#    will be renamed to have the valid date/time of the model state.
#
#-----------------------------------------------------------------------------
# These are the DART namelist variables that must have these values.
# &io_filenames_nml
#   rpointer         = .true.
#   rpointer_file    = 'input_filelist.txt'
# /
# &filter_nml
#   adv_ens_command          = "./advance_model.csh"
#   async                    = 2
#   direct_netcdf_read       = .true.
#   direct_netcdf_write      = .true.
#   init_time_days           = -1
#   init_time_seconds        = -1
#   obs_seq_in_file_name     = "obs_seq.out"
#   start_from_restart       = .true.
# /
# 
#-----------------------------------------------------------------------------
# Python imports here
import os
from netCDF4 import Dataset
from datetime import datetime, timedelta
from namelist_utils import write_namelist, read_namelist



jobname            = 'cm1_filter'
jobid              = os.getpid()

ens_size           = 3

# some systems don't like the -v option to any of the following

copy    = 'cp -pv'
link    = 'ln -sv'
remove  = 'rm -rf'

#-----------------------------------------------------------------------------
# Define the directory that will hold the experiment (i.e. centraldir).
# All the necessary resources get staged here.
# All the run-time input should be defined/modified here.
#-----------------------------------------------------------------------------

centraldir = '{home}/temp/{JOBNAME}/job_{JOBID}'.format(**locals())
centraldir = '/home/disk/jabba/lmadaus/nobackup/osse_ci/{:s}/job_{:d}'.format(jobname, jobid)

print(jobname, jobid, 'centraldir ==', centraldir)

os.system('mkdir -p {:s}'.format(centraldir))
os.chdir(centraldir)

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
# dartdir      The location of the DART cm1 model directory
# cm1dir       The location of the cm1 executable
# ensembledir  The location of the initial ensemble of cm1 files
#-----------------------------------------------------------------------------

dartdir = '/home/disk/pvort/lmadaus/nobackup/cm1/DART_CM1/models/cm1'
cm1dir = '/home/disk/pvort/lmadaus/nobackup/cm1/r18/cm1r18/run'
truthdir = '/home/disk/pvort/lmadaus/nobackup/cm1/DOMAINS/kdvn_ensemble/mems/truth'
ensembledir = '/home/disk/jabba/lmadaus/nobackup/osse_ci/kdvn_ensemble/mems'

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the cm1 executable, control files, and data files.
#-----------------------------------------------------------------------------

os.system('{copy} {dartdir}/work/filter                       . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/advance_time                 . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/input.nml  		          . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/shell_scripts/advance_model.csh   . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/obs_seq.out                  . || exit 1'.format(**locals()))

os.system('{copy} {cm1dir}/cm1.exe              	          . || exit 1'.format(**locals()))

# filter will read the contents of input_filelist.txt for the names
# of the CM1 restart files. After each model advance these same
# filenames will be queried for the new model state.
# So, the actual file associated with this name must be updated.
# advance_model.csh is used by both perfect_model_obs and filter.
# Each model advance is performed in a unique model directory that
# has an instance or ensemble member number in the directory name.

os.system('{remove} input_filelist.txt'.format(**locals()))
os.system('touch input_filelist.txt')

instance = 1
while ( instance <= ens_size ):

    # filter will read the contents of input_filelist.txt for
    # the file containing the most up-to-date model state.
    # Each ensemble member (instance) will get advanced in its own directory.

    filename = 'cm1out_rst_000001.nc'
    dirname = 'dir_model{:03d}'.format(instance)
    os.system('mkdir {:s}'.format(dirname))


    os.system('{copy} {ensembledir}/m{instance}/cm1out_rst_000001.nc {dirname}/{filename} || exit 2'.format(**locals()))
    os.system('{copy} {ensembledir}/m{instance}/input_sounding       {dirname} || exit 2'.format(**locals()))
    os.system('{copy} {ensembledir}/m{instance}/LANDUSE.TBL          {dirname} || exit 2'.format(**locals()))
    os.system('{copy} {ensembledir}/m{instance}/namelist.input       {dirname"/namelist.input.template || exit 2'.format(**locals()))

    # Enforce certain assumptions.
    # 1) we are always starting from a netCDF file ... irst == 1
    # 2) restart always be named cm1_our_rst_000001.nc ... rstnum == 1
    # 3) the run_time value will be a dummy string to facilitate
    #    the use of 'sed' in the advance_model.csh script.
    # Need to modify namelist.input to be a template namelist
    # Read the namelist
    nmld = read_namelist('{:s}/namelist.input.template'.format(dirname))
    nmld['param2']['irst'] = 1
    nmld['param2']['rstnum'] = 1
    # Figure out if this is an environment variable
    cm1_forecast_length = os.environ.get('CM1_FORECAST_LENGTH')
    if cm1_forecast_length not in ['', None]:
        nmld['param1']['run_time'] = cm1_forecast_length
        nmld['param1']['rstfrq'] = cm1_forecast_length

    # Write the new namelist
    write_namelist(nmld, '{:s}/namelist.input'.format(dirname))

    # Might as well put the known ensemble size in the input.nml
    """
    dart_nmld = read_namelist('input.nml')
    dart_nmld['filter_nml']['ens_size'] = ens_size
    dart_nmld['filter_nml']['num_output_obs_members'] = ens_size
    dart_nmld['filter_nml']['num_output_state_members'] = ens_size
    write_namelist(dart_nmld, 'input.nml')
    """

    # append the name of the restart file into the list
    with open('input_filelist.txt','a') as outfile:
        outfile.write('{:s}/{:s}\n'.format(dirname, filename))

    instance += 1

# When filter starts, it needs to read the grid information.
# The DART/CM1 interface specifies the file for the grid information by
# the input.nml:&model_nml:cm1_template_file  variable, which defaults
# to cm1out_rst_000001.nc
#
# Consequently, we can get simply link to any model state.
# Until such time as the 'time' variable is correct, we also need a namelist.input

os.system('{link} dir_model_001/cm1out_rst_000001.nc  cm1out_rst_000001.nc || exit 3'.format(**locals()))
os.system('{copy} dir_model_001/namelist.input.template     namelist.input || exit 3'.format(**locals()))

print("")
print("cd ${centraldir}")
print("Make sure that input.nml contents are correct.")
print("Make sure the obs_seq.in and the model state are appropriate.")
print("Make sure that namelist.input contents are correct.")
print("Make sure that input_filelist.txt contains the right file and the file exists.")
print("Launch (mpirun?) ./filter")
print("")

exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

