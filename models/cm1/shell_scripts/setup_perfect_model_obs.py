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
# Modified for Python (Sept. 2015) Luke Madaus, University of Washington
#-----------------------------------------------------------------------------
#
# This script is meant to provide an example of what is necessary to set up
# a very simple experiment that advances a single instance of CM1 (i.e. the 'truth')
# and harvests a set of synthetic observations from the truth.
#
# ASSUMPTIONS/REQUIREMENTS
#
# *) CM1 is required to be compiled *without* MPI so we can run in async==2
#
# *) CM1 is required to use netCDF restart files.
#
# *) A single CM1 model state will be available as well as a 
#    namelist.input, LANDUSE.TBL, etc.
#
# *) The DART obs_seq.in file will totally define the duration
#    of the truth run.
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
# &perfect_model_obs_nml
#   adv_ens_command          = "./advance_model.csh" # Will overwrite this in script LEM
#   async                    = 2
#   direct_netcdf_read       = .true.
#   direct_netcdf_write      = .true.
#   init_time_days           = -1
#   init_time_seconds        = -1
#   obs_seq_in_file_name     = "obs_seq.in"
#   start_from_restart       = .true.
# /
# 
#-----------------------------------------------------------------------------
from __future__ import print_function, division
import os
from namelist_utils import read_namelist, write_namelist 


jobname    = 'cm1_perfect'
jobid      =  os.getpid()

# Not sure if this variable is set?
cm1_forecast_length =  os.environ.get('CM1_FORECAST_LENGTH')
# Default to one hour
if cm1_forecast_length in ['',None]:
    cm1_forecast_length == 3600.0

# some systems don't like the -v option to any of the following
copy       = 'cp -pv'
link       = 'ln -sv'
remove     = 'rm -rf'

#-----------------------------------------------------------------------------
# Define the directory that will hold the experiment (i.e. CENTRALDIR).
# All the necessary resources get staged here.
# All the run-time input should be defined/modified here.
#-----------------------------------------------------------------------------
home       = os.environ.get('HOME')
#centraldir = os.path.join([home,'temp',jobname,'job_{:d}'.format(jobid)])
centraldir = '/home/disk/jabba/lmadaus/nobackup/osse_ci/cm1_perfect/job_{:d}'.format(jobid)

print(jobname, jobid, "CENTRALDIR ==",  centraldir)

os.system('mkdir -p {:s}'.format(centraldir))
os.chdir(centraldir)

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
# DARTDIR      The location of the DART cm1 model directory
# CM1DIR       The location of the cm1 executable
# TRUTHDIR     The location of the initial model true state
#-----------------------------------------------------------------------------
"""
dartdir = os.path.join(home,'svn','DART','cm1','models','cm1')
cm1dir = os.path.join(home,'svn','DART','cm1','models','cm1','src','cm1r18','run')
truthdir = os.path.join(home,'temp','CM1')
"""
dartdir = '/home/disk/pvort/lmadaus/nobackup/cm1/DART_CM1/models/cm1'
cm1dir = '/home/disk/pvort/lmadaus/nobackup/cm1/r18/cm1r18/run'
truthdir = '/home/disk/pvort/lmadaus/nobackup/cm1/DOMAINS/kdvn_ensemble/mems/truth'


#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the cm1 executable, control files, and data files.
#-----------------------------------------------------------------------------

os.system('{copy} {dartdir}/work/perfect_model_obs          . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/advance_time               . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/input.nml  		        . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/shell_scripts/advance_model.py  . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/obs_seq.in                 . || exit 1'.format(**locals()))

os.system('{copy} {cm1dir}/cm1.exe             	            . || exit 1'.format(**locals()))

# perfect_model will read the contents of input_filelist.txt for the names
# of the CM1 restart file. After each model advance this same
# filename will be queried for the new model state.
# So, the actual file associated with this name must be updated.
# advance_model.csh is used by both perfect_model_obs and filter.
# Each model advance is performed in a unique model directory that
# has an instance or ensemble member number in the directory name.
# So, even though there is only 1 instance for a perfect_model_obs
# experiment, there must still be a directory name compatible with
# the advance_model logic.

filename = "cm1out_rst_000001.nc"
dirname = "dir_model_001"

os.system('mkdir {:s}'.format(dirname))

os.system('{copy} {truthdir}/cm1out_rst_000001.nc {dirname}/{filename} || exit 2'.format(**locals()))
os.system('{copy} {truthdir}/input_sounding       {dirname}            || exit 2'.format(**locals()))
os.system('{copy} {truthdir}/LANDUSE.TBL          {dirname}            || exit 2'.format(**locals()))
os.system('{copy} {truthdir}/namelist.input       {dirname}/namelist.input.template || exit 2'.format(**locals()))

# Enforce certain assumptions.
# 1) we are always starting from a netCDF file ... irst == 1
# 2) restart always be named cm1_out_rst_000001.nc ... rstnum == 1
# 3) the run_time value will be a dummy string
# Python version --use namelist utilities
nmld = read_namelist(os.path.join(dirname, 'namelist.input.template'))
nmld['param2']['irst'] = 1
nmld['param2']['rstnum'] = 1
nmld['param1']['run_time'] = cm1_forecast_length
nmld['param1']['rstfrq'] = cm1_forecast_length
write_namelist(nmld, os.path.join(dirname, 'namelist.input.template'))

os.system('{remove} input_filelist.txt'.format(**locals()))

# Python -- write to input_filelist.txt in with
# statement.  This will automatically close
# the file when we are out of the block
with open('input_filelist.txt','w') as outfile:
    outfile.write(os.path.join(dirname,filename)+'\n')


# When perfect_model_obs starts, it needs to read the grid information.
# The DART/CM1 interface specifies the file for the grid information by
# the input.nml:&model_nml:cm1_template_file  variable, which defaults
# to cm1out_rst_000001.nc
#
# Consequently, we can get simply link to the true model state.

os.system('{link} {dirname}/{filename}                 cm1out_rst_000001.nc || exit 3'.format(**locals()))
os.system('{copy} {dirname}/namelist.input.template    namelist.input       || exit 3'.format(**locals()))

print("")
print("cd {:s}".format(centraldir))
print("Make sure that input.nml contents are correct.")
print("Make sure the obs_seq.in and the model state are appropriate.")
print("Make sure that namelist.input contents are correct.")
print("Make sure that input_filelist.txt contains the right file and the file exists.")
print("Launch ./perfect_model_obs")
print("")

exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

