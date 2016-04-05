#!/bin/csh -f

# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This interactive script builds and submits a series of dependent jobs that runs 
# a multi-cycle DART experiment where multiple model advances and assimilations
# are done in a single batch job, and the 'job dependency' feature of LSF is used 
# to sequence multiple batch jobs. 
#
# Running this way is a bit more complicated than the basic scripts but 
# for a long experiment will be much cheaper and should have faster throughput.
# See the README file for more details on how to use these scripts.
#
# This script is intended to work only with the LSF batch system.  It would
# need to be heavily adapted for other batch systems.

# (to actually use this script, remove the following lines once you
# have reviewed how these scripts work.)
echo "HEY!  This is a script intended only for advanced users "
echo "      on NCAR's yellowstone computer."
echo "      Please consult us before using it, since it does not have the level of "
echo "      error checking and commenting which other DART scripts have. "
exit 218                                                                            
 
# ==============================================================================

# Run this from the CESM $CASE directory

# Each 'job' runs RESUBMIT+1 forecasts+assimilation cycles.
#    This is defined for CESM by the RESUBMIT variable in env_run.xml,
#    and is set in this script, and then reset by each job at its beginning.
set RESUBMIT = 9

# 'jobs_loop' jobs can be run before the $RUNDIR disk fills up
#    and the output must be archived/purged.  
set jobs_loop = 1

# Long assimilations may require several archivings.  
#    This is defined by archive_loop.
set archive_loop = 1

@ cycles = ($RESUBMIT + 1) * $jobs_loop * $archive_loop
echo "The total number of forecasts submitted by this script is $cycles"

# Set RESUBMIT for the first job, to make $CASE.run_cycles start with the 
# right number of forecasts.  It will take care of itself after this.
./xmlchange RESUBMIT=$RESUBMIT
./xmlchange CONTINUE_RUN=TRUE

source ./Tools/ccsm_getenv  || exit 14
echo RESUBMIT is now $RESUBMIT

# First job has no dependency built in.
setenv BATCHSUBMIT_DEP `echo $BATCHSUBMIT`

# Loop over the number of archivings which will be required.
set a_l = 1
while ($a_l <= $archive_loop)
   echo "archive loop $a_l"

   # Loop over the number of jobs which fill up the scratch space,
   # and necessitate archiving.
   # Each job will start with the RESUBMIT value set in this script.
   set j_l = 1
   while ($j_l <= $jobs_loop)
      echo "   job loop $j_l"
      echo "      ${BATCHSUBMIT_DEP} ./${CASE}.run_cycles" 

      echo "${BATCHSUBMIT_DEP} ./${CASE}.run_cycles" >! templar
      source templar >! out

      set i = `grep Job out`
      set jobid = `echo $i | cut -d'<' -f2 | cut -d'>' -f1`

      # The next job will start if the previous job state is DONE (not EXIT).
      setenv BATCHSUBMIT_DEP  `echo 'bsub -w "done(' $jobid ')" <'`
      @ j_l++
   end

   # ONLY submit the post run cleanup, archive and long_term archive
   # if completion of the last job (all of its RESUBMITs) is successful.

   setenv BATCHSUBMIT_DEP  `echo 'bsub -w "done(' $jobid ')" <'`
   echo "   ${BATCHSUBMIT_DEP} ./archive_cycles.csh"  

   echo "${BATCHSUBMIT_DEP} ./archive_cycles.csh"  >! templar
   source templar >! out

   set i = `grep Job out`
   set jobid = `echo $i | cut -d'<' -f2 | cut -d'>' -f1`

   # submit the next $CASE.run_cycles only if the current archive_cycles.csh is "done"
   setenv BATCHSUBMIT_DEP  `echo 'bsub -w "done(' $jobid ')" <'`

   @ a_l++
end

if ($status == 0) then
   rm templar out
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
