#!/bin/csh
#
# This is a script to create a series of conditional jobs.
# When one job terminates normally, the next one enters the queue.
# The LSF queueing system does most of the work, this script
# simply builds the appropriate LSF headers for the series of jobs.
#
# This script is meant to be run interactively.
# It will cause a series of jobs to get submitted to the queue.
#
# TJH -- 2 June 2006 
#----------------------------------------------------------------------

# Create a unique identifier for this experiment,
# Create a unique counter to facilitate the creation of sequential,
# unique job names -- needed by LSF.

set experiment = smoother
set jobnumber = 1
set njobs = 10

#----------------------------------------------------------------------
# Loop around the number of obs sequence files?
#----------------------------------------------------------------------

while ( $jobnumber <= $njobs ) 

   set jobname = $experiment$jobnumber

   #-------------------------------------------------------------------
   # Build the batch system directives, put into a file.
   #-------------------------------------------------------------------

   rm -f ${jobname}.lsf

   echo "#!/bin/csh"                > ${jobname}.lsf
   echo "#BSUB -J $jobname"        >> ${jobname}.lsf
   echo "#BSUB -o $jobname.%J.log" >> ${jobname}.lsf
   echo "#BSUB -q regular"         >> ${jobname}.lsf
   echo "#BSUB -n 1"               >> ${jobname}.lsf

   if ($jobnumber > 1) then
      @ previousjobnumber = $jobnumber - 1
      set previousjobname = $experiment$previousjobnumber
      echo "#BSUB -w done($previousjobname)" >> ${jobname}.lsf
   endif

   #-------------------------------------------------------------------
   # append the executable statements to the file, now it is a script.
   #-------------------------------------------------------------------

   copy input and restart files to the directory

   echo "cd /image/home/khare/DART/models/lorenz_96/work" >> ${jobname}.lsf
   echo "mkdir $jobname"                                  >> ${jobname}.lsf
   echo "./perfect_model_obs"                             >> ${jobname}.lsf
   echo "./smoother"                                      >> ${jobname}.lsf
   echo "cp -pv obs_seq.final $jobname"                   >> ${jobname}.lsf

   #-------------------------------------------------------------------
   # execute the script
   #-------------------------------------------------------------------

   echo " "
   echo "Script $jobnumber looks like: "
   cat ${jobname}.lsf
   echo " "

#   bsub < ${jobname}.lsf

   @ jobnumber++

end
