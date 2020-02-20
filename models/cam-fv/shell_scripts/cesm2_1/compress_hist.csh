#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Compress just the history files from all of the components,
# one type (h0,...) at a time, which limits the PEs needed
# to the size of the ensemble.

# ------------------------------------------------------------------------------

if ($#argv != 5) then
   echo "Usage: In the directory containing the files to be processed:"
   echo "   call with exactly 5 arguments or submit as a batch job with 0 arguments:"
   echo '   ${scr_dir}/compress_hist.csh command YYYY-MM-DD-SSSS "sets" "types" "stages"'
   echo '   where '
   echo '   sets   = 1 or more of {clm2 cpl cam cice hist dart} to compress, separated by spaces'
   echo '   types  = 1 or more CESM+DART "file types" (h0, hr2x,...) to compress'
   echo '   stages = 1 or more of stages {input, forecast, postassim, output} to compress.'
   echo ' -OR-'
   exit 17
endif

# Environment variables ($data_*) from the calling script should be available here.

set comp_cmd      = $1
set ymds          = $2
set sets          = ($3)
set types         = ($4)
set stages        = ($5)

set cmd = `echo $comp_cmd | cut -d' ' -f1`
if ($cmd == 'gzip') then
   set ext = ''
else if ($cmd == 'gunzip') then
   set ext = '.gz'
else
   echo "ERROR: unrecognized command $cmd.  Don't know which extension to use"
   exit 27
endif

echo "In compress_hist.csh:"
echo "   (arg1) comp_cmd      = $comp_cmd"
echo "   (arg2) date          = $ymds"
echo "   (arg3) sets          = $sets"
echo "   (arg4) types         = $types"
echo "   (arg5) stages        = $stages"
echo "   data_CASE     = $data_CASE"
echo "   data_CASEROOT = $data_CASEROOT"
echo "   data_NINST    = $data_NINST"
# echo "   data dir      = $data_dir"

# ------------------------------------------------------------------------------
# Fail if there are leftover error logs from previous compression.csh executions.

ls *.eo > /dev/null
if ($status == 0) then
   echo "ERROR; Existing compression log files: *.eo.  Exiting"
   exit 37
endif

# --------------------------
# Environment and commands.


setenv MPI_SHEPHERD true

setenv date 'date --rfc-3339=ns'

# ------------------------------------------------------------------------------
# Create the command file where each line is a separate command, task, operation, ....

foreach type ( $types )

   \rm -f mycmdfile
   touch mycmdfile

   # 'task' is incremented continuously over all files; components, members, etc.
   # 'task' is a running counter of jobs in mycmdfile.
   set task = 0
   
   foreach comp ( $sets )
   echo "comp = $comp"
   switch ($comp)
      case {clm2,cpl,cam,cice,mosart}:
         set i=1
         while ( $i <= $data_NINST)
            # E.g. CAM6_80mem.cice_0001.r.2010-07-15-00000.nc
            set file_name = \
                `printf "%04d/%s.%s_%04d.%s.%s.nc%s" $i $data_CASE $comp $i ${type} $ymds $ext`
            # If the expected file exists, add the compression command 
            if (-f $file_name) then
               @ task++
               echo "$comp_cmd $file_name &> compress_${task}.eo " >> mycmdfile
            else if (-f ${file_name}.gz) then
               # Kluge to get around situations where an earlier job compressed the file,
               # but failed for some other reason, so it's being re-run.
               echo "$file_name already compressed"
            else
               echo 'ERROR: Could not find "'$file_name'" to compress.'
               exit 47
            endif
   
            @ i++
         end
      breaksw

      default:
         breaksw
   endsw
   end

   # ------------------------------------------------------------------------------

   echo "Before launching mycmdfile"

   $date 
   
   # CHECKME ... make sure $task is <= the number of MPI tasks in this job.
   
   if ($task > 0) then
      if ($?PBS_O_WORKDIR) then
         mpiexec_mpt -n $task ${data_CASEROOT}/launch_cf.sh ./mycmdfile
         set mpi_status = $status
      else if ($?SLURM_SUBMIT_DIR) then
         mpirun      -n $task ${data_CASEROOT}/launch_cf.sh ./mycmdfile
         set mpi_status = $status
      endif
   
      echo "mpi_status = $mpi_status"
   else
      echo "No compression to do"
      exit 0
   endif
   
   # Check the statuses?
   if ( -f compress_1.eo ) then
      grep $cmd *.eo
      # grep failure = compression success = "not 0"
      set gr_stat = $status
   #    echo "gr_stat when eo file exists = $gr_stat"
   else
      # No eo files = failure of something besides g(un)zip.
      set gr_stat = 0
   #    echo "gr_stat when eo file does not exist = $gr_stat"
   endif
   
   if ($gr_stat == 0) then
      echo "compression failed.  See .eo files with non-0 sizes"
      echo "Remove .eo files after failure is resolved."
      exit 197
   else
      # Compression worked; clean up the .eo files and mycmdfile
      \rm -f *.eo  mycmdfile
   endif
   
# types
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
