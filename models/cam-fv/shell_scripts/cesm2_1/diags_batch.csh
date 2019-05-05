#!/bin/tcsh

#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -W 1:30
#BSUB -q caldera
#BSUB -P P86850054
#BSUB -N
#BSUB -u raeder@ucar.edu
#BSUB -o diags_batch.%J
#BSUB -e diags_batch.%J
#BSUB -J diags_batch
#-----------------------------------------
#PBS  -N diags_batch
#PBS  -A P86850054
#PBS  -q share
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
#PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=02:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M raeder@ucar.edu
# Send standard output and error to this file.
# It's helpful to use the $casename here.
#PBS  -o diags_batch.eo
#PBS  -j oe 
#--------------------------------------------

#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source: /home/thoar/CVS.REPOS/DART/models/cam/shell_scripts/diags.csh,v $
# $Name:  $

# DART source directory on this machine

# set echo verbose

# obs_diag can now take a namelist argument that specifies a file
# containing a list of obs_seq filenames.
if ($#argv == 0) then
   if ($?LS_SUBCWD) then
      cd $LS_SUBCWD
   else if ($?PBS_O_WORKDIR) then
      cd $PBS_O_WORKDIR
   endif
   set case = something_meaningful
   set diag_dir = Diags.${case}_NTrS_2017.1.1-2017.1.8H0_s0
   set obs_dir = ../Obs_${case}_noAIRS
else
   set case     = $1
   set obs_dir  = ../$2
   set diag_dir = $3
   
endif

# These things should be gathered from env_*.xml files in CASEROOT.
set DART = ~/DART/reanalysis
set cam = 'cam-fv'
set machine = '_casper'
# Cheyenne;  set machine = ''

# Use big endian obs_diag for POP_force output from IBM
# set endian = '_big_endian'
set endian = ' '

if (! -d $obs_dir:t) then
   echo "Missing obs_dir"
   exit 10
# else
# Done in calling script    ln -s ${obs_dir}.list obs.list
endif

if (! -d $diag_dir) then
   mkdir $diag_dir
   cd $diag_dir
else
   echo "$diag_dir exists; choose another name"
   exit 20
endif

# set echo verbose
pwd

# Done in obs_seq_tool_series.csh
# ls ${obs_dir}/*obs_seq_final*common_3_AIRS >! obs.list 
# Replacement
ln -s ${obs_dir}.list obs.list

cp ../input.nml .

echo "Running ${DART}/models/${cam}/work/obs_diag${endian}"
${DART}/models/${cam}/work${machine}/obs_diag${endian} >&! obs_diag.out 

exit
