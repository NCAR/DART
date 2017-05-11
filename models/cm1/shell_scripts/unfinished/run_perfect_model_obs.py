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
# Top level script to generate observations and a TRUE state.
#
# This script is designed to be submitted as a batch job but may be run from 
# the command line (as a single thread) to check for file motion, etc.
# If running interactively, please comment out the part that actually runs filter.
#
#-----------------------------------------------------------------------------
#
#BSUB -J cm1_perfect
#BSUB -o cm1_perfect.%J.log
#BSUB -P P3507xxxx
#BSUB -q premium
#BSUB -n 1
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu
from __future__ import print_function, division
import os
from datetime import datetime, timedelta
from netCDF4 import Dataset
from namelist_utils import read_namelist, write_namelist


#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if (os.environ.get('LSB_HOSTS') not in ['', None]):

    originaldir = os.environ.get('LS_SUBCWD')
    jobname     = os.environ.get('LSB_JOBNAME')
    jobid       = os.environ.get('LSB_JOBID')
    myqueue     = os.environ.get('LSB_QUEUE')
    myhost      = os.environ.get('LSB_SUB_HOST')

else:

    #-------------------------------------------------------------------
    # You can run this interactively to check syntax, file motion, etc.
    #-------------------------------------------------------------------

    originaldir = os.getcwd()
    jobname     = 'cm1_perfect'
    jobid       = os.getpid()
    myqueue     = 'Interactive'
    myhost      = os.environ.get('HOST')



#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------
print("")
print(jobname, ' (', jobid, ') submitted   from ', originaldir)
print(jobname, ' (', jobid, ') submitted   from ', myhost)
print(jobname, ' (', jobid, ') running in queue ', myqueue)
print(jobname, ' (', jobid, ') running       on ', myhost)
print(jobname, ' (', jobid, ') started   at ',  datetime.now())
print("")

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#----------------------------------------------------------------------
centraldir = '/home/disk/jabba/lmadaus/nobackup/osse_ci/{:s}/job_{:d}'.format(jobname,jobid)

os.system('mkdir -p {:s}'.format(centraldir))
os.chdir(centraldir)

myname = 'run_perfect_model_obs'

# some systems don't like the -v option to any of the following 

remove = 'rm -rvf'
copy   = 'cp -v'
move   = 'mv -fv'
link   = 'ln -s'

print(jobname, jobid, 'centraldir == ', centraldir)

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
# DARTDIR      The location of the DART cm1 model directory
# cm1DIR    The location of the cm1 executable
# ENSEMBLEDIR  The location of the initial ensemble of cm1 files
# EXPERIMENT   The (safe) location for the results of this run.
#-----------------------------------------------------------------------------
"""
set     DARTDIR = /home/disk/pvort/nobackup/${USER}/cm1/DART_CM1/models/cm1
set   cm1DIR = /home/disk/pvort/nobackup/${USER}/cm1/r18/cm1r18/run
set ENSEMBLEDIR = /home/disk/pvort/nobackup/${USER}/cm1/DOMAINS/kdvn_ensemble/mems/truth
set  EXPERIMENT = /home/disk/pvort/nobackup/${USER}/cm1/DART_CM1/models/cm1/work
"""
"""
dartdir = os.path.join(home,'svn','DART','cm1','models','cm1')
cm1dir = os.path.join(home,'svn','DART','cm1','models','cm1','src','cm1r18','run')
ensembledir = os.path.join(home,'temp','CM1')
"""
dartdir = '/home/disk/pvort/lmadaus/nobackup/cm1/DART_CM1/models/cm1'
cm1dir = '/home/disk/pvort/lmadaus/nobackup/cm1/r18/cm1r18/run'
ensembledir = '/home/disk/pvort/lmadaus/nobackup/cm1/DOMAINS/kdvn_ensemble/mems/truth'
experiment = '/home/disk/pvort/lmadaus/nobackup/cm1/DOMAINS/kdvn_ensemble/mems/truth'

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the cm1 executable, control files, and data files.
# The cm1 initial conditions are in the next block.
#-----------------------------------------------------------------------------

os.system('{copy} {dartdir}/work/perfect_model_obs          . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/advance_time               . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/input.nml                  . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/shell_scripts/advance_model.py  . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/shell_scripts/namelist_utils.py . || exit 1'.format(**locals()))
os.system('{copy} {dartdir}/work/obs_seq.in                 . || exit 1'.format(**locals()))

os.system('{copy} {cm1dir}/cm1.exe                          . || exit 1'.format(**locals()))


os.system('{copy} {ensembledir}/cm1out_rst_000001.nc        . || exit 2'.format(**locals()))
os.system('{copy} {ensembledir}/input_sounding              . || exit 2'.format(**locals()))
os.system('{copy} {ensembledir}/LANDUSE.TBL                 . || exit 2'.format(**locals()))
os.system('{copy} {ensembledir}/namelist.input              . || exit 2'.format(**locals()))



# Enforce certain assumptions.
# 1) we are always starting from a netCDF file ... irst == 1
# 2) restart always be named cm1_out_rst_000001.nc ... rstnum == 1
# 3) the run_time value will be a dummy string
# Python version --use namelist utilities
# read the namelist
nmld = read_namelist('namelist.input')

# Check for environment variable that
# contains the current forecast length
# Otherwise leave as is
if os.environ.get('CM1_FORECAST_LENGTH') not in [None, '']:
    cm1_forecast_length = float(os.envion.get('CM1_FORECAST_LENGTH'))
    nmld['param1']['run_time'] = cm1_forecast_length
    nmld['param1']['rstfrq'] = cm1_forecast_length
                              
# Set the other namelist requirements and write                               
nmld['param2']['irst'] = 1
nmld['param2']['rstnum'] = 1
write_namelist(nmld, 'namelist.input')



# Python -- write to input_filelist.txt in with
# statement.  This will automatically close
# the file when we are out of the block
# Only have one ens member to worry about here, so link directly
os.system('{remove} input_filelist.txt'.format(**locals()))
with open('input_filelist.txt','w') as outfile:
    outfile.write('dir_model_001/cm1out_rst_000001.nc\n')

# While we're here, make a copy of the restart file for
# advance model with the time of the current restart
# file attached to it
# Read the simulation start (base time) from the namelist
nmldtime = nmld['param11']
startdt = datetime(nmldtime['year'], nmldtime['month'], nmldtime['day'],\
                   nmldtime['hour'], nmldtime['minute'], nmldtime['second'])
# Read the current time of the netcdf restart file
# This opens a pointer to the file
with Dataset('cm1out_rst_000001.nc','r') as dset:
    # Read the time variable for the duration
    # of the with-block
    cur_time = int(dset.variables['time'][0])
# Add this to the startdt for the currentdt
currentdt = startdt + timedelta(seconds=cur_time)

# Advance model needs instance number in most cases, add one here
os.system('{:s} cm1out_rst_000001.nc cm1out.001.{:%Y%m%d%H%M%S}.nc'.format(copy, currentdt))

# And for our pointer file, need to have this file linked into a member directory
# At the beinning.  Only one ens member, so easy to create:
if not os.path.exists('dir_model_001'):
    os.system('mkdir dir_model_001')
os.chdir('dir_model_001')
os.system('{:s} ../cm1out_rst_000001.nc .'.format(link))
os.chdir(centraldir)

#-----------------------------------------------------------------------------
# Run perfect_model_obs ... harvest the observations to populate obs_seq.out
# model_mod expects a generic name // advance_model.csh expects a filename
# with the ensemble member ID tacked on - must provide both.
#-----------------------------------------------------------------------------

os.system('./perfect_model_obs || exit 3')

#-----------------------------------------------------------------------------
# At this point, all the restart,diagnostic files are in the run/CENTRALDIR.
# You may want to move them to someplace more 'permanent'.
#
# TJH: At this point, the output files have pretty 'generic' names.
# The files could be archived with the assimilation date in their name.
#-----------------------------------------------------------------------------

# os.system('{move} cm1_s.nc.0001           {experiment}/perfect/cm1_s.nc'.format(**locals()))
# os.system('{move} cm1_restart_p.nc.0001   {experiment}/perfect/cm1_restart_p.nc'.format(**locals()))
# os.system('{move} cm1.nml                 {experiment}/perfect'.format(**locals()))
# os.system('{move} obs_seq.out                {experiment}/perfect'.format(**locals()))
# os.system('{move} True_State.nc              {experiment}/perfect'.format(**locals()))

# os.system('{move} cm1_out_1               {experiment}/perfect/cm1_out'.format(**locals()))
# os.system('{move} dart_log.out               {experiment}/perfect'.format(**locals()))
# os.system('{move} dart_log.nml               {experiment}/perfect'.format(**locals()))
# Good style dictates that you save the scripts so you can see what worked.

# os.system('{copy} input.nml                  {experiment}/DART'.format(**locals()))
# os.system('{copy} *.csh                      {experiment}/DART'.format(**locals()))
# os.system('{copy} $myname                    {experiment}/DART'.format(**locals()))

print("{jobname} ({jobid}) finished at".format(**locals()), datetime.now())
print("These are the files in the run directory at completion:")
os.system('ls -lrt')

exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

