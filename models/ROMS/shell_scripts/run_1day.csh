#!/bin/tcsh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#=============================================================================
#
#
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluefire'
#
# the normal way to submit to the queue is:    bsub < roms_ensemble.csh
#
# an explanation of the most common directives follows:
# -J Job_name
# -o STDOUT_filename
# -e STDERR_filename
# -P account_code_number
# -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
# -n number of MPI processes (not nodes)
# -W hh:mm  wallclock time (required on some systems)
#=============================================================================
#BSUB -J roms_ensemble
#BSUB -o roms_ensemble.%J.log
#BSUB -P P86850054
#BSUB -q small
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu
#
#=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
## PBS is used on the CGD Linux cluster 'bangkok'
## PBS is used on the CGD Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub run_roms_ensemble.csh
##
## an explanation of the most common directives follows:
## -N     Job name
## -r n   Declare job non-rerunable
## -e <arg>  filename for standard error
## -o <arg>  filename for standard out
## -q <arg>   Queue name (small, medium, long, verylong)
## -l nodes=xx:ppn=2   requests BOTH processors on the node. On both bangkok
##                     and calgary, there is no way to 'share' the processors
##                     on the node with another job, so you might as well use
##                     them both. (ppn == Processors Per Node)
#=============================================================================
#PBS -N roms_ensemble
#PBS -r n
#PBS -e roms_ensemble.err
#PBS -o roms_ensemble.log
#PBS -q dedicated
#PBS -l nodes=10:ppn=2

#>todo FIXME : make this automagic for the setup script

cd /glade/scratch/thoar/romstest2

@ instance = 0
foreach INSTANCE_DIRECTORY ( instance_???? )

   @ instance++

   cd ${INSTANCE_DIRECTORY}

   rm -f log_$instance.txt

   mpirun.lsf ../oceanM ocean.in |& tee log_$instance.txt

   # Check for successful completion - log file should have something like:
   # ROMS/TOMS: DONE... Tuesday - April 26, 2016 -  7:34:13 PM
   grep "ROMS/TOMS: DONE" log_$instance.txt > /dev/null
   if ($status != 0) then
      echo "ROMS instance $instance FAILED."
      echo "ROMS instance $instance FAILED."
      echo "ROMS instance $instance FAILED."
      exit 1
   endif

   # tag the output with the model time 
   # Create a ROMS_POSTERIOR file that will be updated by DART.
   # Save off a copy of the forecast as the PRIOR.
   set OCEAN_TIME_STRING = `ncdump -v ocean_time roms_dai.nc | grep '^ ocean_time = '`
   set OCEAN_TIME = `echo $OCEAN_TIME_STRING | sed -e "s#[=;a-z_ ]##g"`

   set ROMS_PRIOR     = `printf roms_dai_%04d_%d.nc $instance $OCEAN_TIME`
   set ROMS_POSTERIOR = `printf roms_pos_%04d_%d.nc $instance $OCEAN_TIME`
   set ROMS_OBSFILE   = `printf roms_mod_%04d_%d.nc $instance $OCEAN_TIME`

   \cp -v roms_dai.nc      ${ROMS_POSTERIOR}
   \mv -v roms_dai.nc      ${ROMS_PRIOR}
   \mv -v roms_obs_mod.nc  ${ROMS_OBSFILE}

   echo
   echo "#---------------------------------------------------------------------"
   echo "# ROMS instance $instance completed at "`date`
   echo "#---------------------------------------------------------------------"
   echo

   cd ..

end

#==============================================================================
# Then we run DART on the ensemble of new states
#==============================================================================

# Remove the last set of DART run-time logs - if they exist.
\rm -f dart_log.out dart_log.nml

# Because convert_roms_obs and filter need bits from the ROMS model_mod,
# a (single) ROMS input file is required to satisfy 'static_init_model()' 
# any one will do

ln -sf instance_0001/roms_dai_????_${OCEAN_TIME}.nc roms_input.nc

# Collect all the ROMS_OBSFILEs into a list of input files 
# and then convert them to a single DART observation sequence file.

ls -1 instance_*/roms_mod*_${OCEAN_TIME}.nc  >! precomputed_files.txt

./convert_roms_obs  || exit 2

# 2) collect all the ROMS_RESTARTs into a list of input files for filter
# The io module will error out if the file_list.txt is too short.
# (make sure all instances of ROMS advanced successfully)
# DART (filter) will modify these files in-place.

ls -1 instance_*/roms_pos_????_${OCEAN_TIME}.nc  >! restart_files.txt

mpirun.lsf ./filter || exit 3

#==============================================================================
# Prepare for the next model advance
#==============================================================================

# 1) filter will write out a 'new_dstart.txt' file with the new DSTART
#    that must be inserted into the CENTRALDIR/ocean_wc13.in

\mv -v PriorDiag_sd.nc    Prior_sd.${OCEAN_TIME}.nc
\mv -v PriorDiag_mean.nc  Prior_mean.${OCEAN_TIME}.nc
\mv -v sd.nc              Posterior_sd.${OCEAN_TIME}.nc
\mv -v mean.nc            Posterior_mean.${OCEAN_TIME}.nc
\mv -v obs_seq.final      obs_seq.final.${OCEAN_TIME}

# what inflation file 

# the input.nml:&filter_nml:restart_out_file_name specifies the base filename
# for the updated (the posterior) model state. This base gets appended with
# a 4-digit instance number and then ".nc" This file must be pushed back into
# the appropriate directory. 

@ instance = 0
foreach INSTANCE_DIRECTORY ( instance_???? )
   @ instance++

   cd ${INSTANCE_DIRECTORY}

   set posterior = `head -n $instance ../restart_files.txt | tail -n 1`
   set  filename = $posterior:t

   # use new state as starting point for next advance.
   # We want to preserve the unique posterior but want ROMS
   # to read from the (single) filename in ocean.in INIFILE

   ln -sf ${filename}  roms_input.nc || exit 4

   cd ..

end


#==============================================================================
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

