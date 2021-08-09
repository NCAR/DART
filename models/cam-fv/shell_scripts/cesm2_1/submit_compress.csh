#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# $Id$

#==========================================================================
#
#SBATCH --job-name=compress
#SBATCH -o %x_%j.errout 
#SBATCH -e %x_%j.errout 
# 80 members
# restarts 
#SBATCH --ntasks=320 
# forcing files: #SBATCH --ntasks=405 
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
#SBATCH --account=NCIS0006
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------
#PBS  -N compress.csh
#PBS  -A YOUR_ACCOUNT
#PBS  -q premium
# For restarts:
# #PBS  -l select=9:ncpus=36:mpiprocs=36
# For hist: 5 * 80         = 400  / 36 = 12
# For dart: 1 + 2*(2 + 80) = 165  
#                            645 / 36 = 18
# For rest: 4 * 80         = 320 / 36 =  9
#PBS  -l select=18:ncpus=36:mpiprocs=36
#PBS  -l walltime=00:20:00
#PBS  -o compress.out
#PBS  -j oe 
# Submit from CASEROOT instead of the rest directory,
# so that the job output file doesn't clutter the data directory
# and xml file variables are available.
set DOUT_S_ROOT = `./xmlquery DOUT_S_ROOT --value`

set comp_cmd      = 'gzip '
# set comp_cmd      = 'gzip -k'
set ymds          = 2019-07-01-00000
set data_dir      = ${DOUT_S_ROOT}/rest/${ymds}

# set sets          = (cpl)
# Restarts
set sets          = (clm2 cpl cam cice)
# set types         = ( ha2x1d hr2x ha2x3h ha2x1h ha2x1hi )
set stages        = (none)

cd ${data_dir}

${data_CASEROOT}/compress.csh $comp_cmd $ymds "$sets" "$stages"
