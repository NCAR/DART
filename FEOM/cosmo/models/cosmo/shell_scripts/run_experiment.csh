#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Top level script to run a COSMO assimilation experiment.
#

#=============================================================================
#  initial setup
#=============================================================================

set expname  = bob        # name of this experiment
set ncycles  = 10         # how many assimilation/model advance loops to do
set num_ens  = 32         # how many ensemble members?
set proj_num = xxxxxxxx   # project to charge this job to
set nprocs   = 32         # number of mpi tasks for filter
set num_nodes = 4         # number of nodes for each COSMO model run

set nhours       = 3
set start_time   = 2009-10-01_00:00:00  # oct 1st, 2009, 0Z

set adv_time     = +${nhours}h                        # N hour model advances
set elapsed_time = `printf 00%02d0000 $nhours`        # N hours in DDhhmmss format

set COSMO_DIR = cosmo
set obs_base = /glade/proj3/image/Observations/ACARS/   # where to find obs
set run_dir  = .

set usingsec = 'false'    # sampling error correction option
set usinginf = 'true'     # using adaptive inflation

#set submit   = 'bsub < '   # or qsub, or qsub <, or mpirun -n X
set submit   = 'echo '   # or qsub, or qsub <, or mpirun -n X
set startmpi = mpirun.lsf  # or mpirun or  ...

set reqfiles = "input.nml.template advance_time filter fill_inflation_restart cosmo_to_dart dart_to_cosmo"

foreach i ( $reqfiles )
   if ( ! -e $i ) then
      echo required file $i not found, exiting.
      exit -2
   endif
end

# if you are using sampling error correction, you need a file that
# depends on the number of ensemble members.
if ( $usingsec == 'true' ) then
   if ( ! -e final_full.${num_ens} ) then
      echo sampling error correction file not found: final_full.${num_ens}
      echo get from DART/system_simulation/final_full_precomputed_tables
      echo or generate your own following the instructions in the html page
      echo -3
   endif
endif

# if you are using adaptive inflation and this is the first job step,
# use 'fill_inflation_restart' to initialize a set of restart files.

# make a dummy input.nml if there isn't one, so we can run advance_time
if ( ! -e input.nml ) then
   echo '&utilities_nml' >! input.nml
   echo ' /'             >> input.nml
endif
set currtime = `echo $start_time 0 | ./advance_time`

# set up the common header things that will be used in all job submit scripts
echo "s/JOB_NAME/$expname/g"     >! script1.sed
echo "s/PROJ_NUMBER/$proj_num/g" >> script1.sed
echo "s/NPROCS/$nprocs/g"        >> script1.sed

# if needed, sed anything that changes and do:
#  sed -f script1.sed input.nml.template > input.nml
cp input.nml.template input.nml

#
#=============================================================================
#  main cycle loop
#=============================================================================

set icyc = 1
while ( $icyc < $ncycles )

   echo Starting cycle $icyc of $ncycles
   echo current time: $currtime

   # things that may change with cycles, or just make more sense here
   echo "s;OBSBASE;$obs_base;g"      >! script2.sed
   echo "s/YYYYMMDDhh/$currtime/g"   >> script2.sed
   echo "s/MPIRUN/$startmpi/g"       >> script2.sed
   echo "s/CURRTIME/$currtime/g"     >> script2.sed
   echo "s/STARTTIME/$start_time/g"  >> script2.sed

   cat script1.sed script2.sed >! script.sed

   #==========================================================================
   #  1. convert model output files to filter input files
   #==========================================================================
   
   # args are:  
   #   dir where cosmo files are found
   #   dir where dart files should go
   #   name of cosmo input files
   #   base name of dart output files
   #   first ensemble number
   #   last ensemble number
   #   number of digits for cosmo subdirectory names
   #   number of digits for dart filenames
   ./cosmo_to_dart.sh $COSMO_DIR $rundir lfff${elapsed_time} filter_ics. 1  ${num_ens} 3 4
   
   #==========================================================================
   #  2. set up and run filter to assimilate next set of observations
   #==========================================================================
   
   # collect files you have from a previous run:
   # restart, inflation, etc

   # need to link obs to current time
   ln -s $obs_base/obs_seq${currtime} obs_seq.out
   eval '$submit ./filter'
   
   #==========================================================================
   #  3. convert filter output to model input files
   #==========================================================================
   
   # args are:  
   #   dir where dart files should go
   #   dir where cosmo files are found
   #   base name of dart output files
   #   name of cosmo input files
   #   first ensemble number
   #   last ensemble number
   #   number of digits for dart filenames
   #   number of digits for cosmo subdirectory names
   ./dart_to_cosmo.sh $rundir $COSMO_DIR filter_restart. laf${currtime} 1 $num_ens 4 3
   
   #==========================================================================
   #  4. run N copies of the model here
   #==========================================================================
   
   # args are:
   #  cosmo run directory
   #  current time
   #  number of hours to run
   #  start ensemble number
   #  end ensemble number
   #  number of nodes for each COSMO model run
   ./run_cosmo_ensemble.sh $COSMO_DIR $currtime $nhours 1 $num_ens $num_nodes
   
   #==========================================================================
   #  5. (optional) archive files here
   #==========================================================================
   
   echo add hsi or mss commands here

   # increment time and loop again
   set currtime = `echo $currtime $adv_time | ./advance_time`
   @ icyc ++
end

#=============================================================================
#  bottom of main cycle loop
#=============================================================================

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


