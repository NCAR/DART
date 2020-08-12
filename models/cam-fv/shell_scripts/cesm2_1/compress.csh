#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#PBS  -N compress.csh
#PBS  -A P86850054
#PBS  -q premium
# For restarts:
# #PBS  -l select=9:ncpus=36:mpiprocs=36
# For hist: 6 * 80         = 480
# For dart: 1 + 2*(2 + 80) = 165  
#                            645 / 36 = 18
# For rest: 4 * 80         = 320 / 36 =  9
#PBS  -l select=18:ncpus=36:mpiprocs=36
#PBS  -l walltime=00:20:00
#PBS  -o compress.out
#PBS  -j oe 

# ------------------------------------------------------------------------------
# Purpose:
#
# compresses or uncompresses sets of files from a forecast or assimilation.
#
# ------------------------------------------------------------------------------
#
# Method:
#
# When called from a script (normally assimilate.csh), it compresses the files.
# When submitted as a batch job, it can also uncompress sets of files -  
# IF this script has been configured to use the right metadata to
# construct the expected data directories.
#
# The strategy is to create a cmdfile that contains a separate task on each line
# and dispatch that cmdfile to perform N simultaneous operations. That cmdfile
# has a syntax ( &> ) to put stderr and stdout in a single file.
# PBS requires 'setenv MPI_SHEPHERD true' for the cmdfile to work correctly.
#
# Compression method can depend on the file type, and in the future may include 
# lossy compression. This script is most often called by assimilate.csh, but can 
# be run as a batch job.
#
# Assimilate.csh runs this in 2 places:
# 1) Every cycle: 
#    +  all the cpl history (forcing) files.
#    +  DART output
#       >  stages of state files  
#             mean, sd  (no instance number)
#       >  obs_seq.final (no instance number)
#       >  Note: inflation files are never compressed.
# 2) Before archiving a restart set to archive/rest; all large restart files.
# ------------------------------------------------------------------------------

if ($#argv == 5) then
   # Called from assimilate.csh (or other script).
   set comp_cmd      = 'gzip'
   set case_name     = $1
   set ymds          = $2
   set ensemble_size = $3
   set sets          = ($4)
   set stages        = ($5)
   set data_dir      = '.'

else if ($#argv == 0) then
   # Edit these and run as a batch job.
   # 'sets' performs better when ordered by decreasing size (clm2 cpl cam cice hist dart)
   set comp_cmd      = 'gunzip'
   set case_name     = CESM2_1_80_3node
   set ymds          = 2010-07-17-64800
   set ensemble_size = 80
   set sets          = (hist dart)
   # set sets          = (clm2 cpl cam cice)
   set stages        = (preassim output)
   # set data_dir      = /glade/scratch/${USER}/${case_name}/archive/rest/${ymds}
   set data_dir      = /glade/scratch/${USER}/${case_name}/run

else 
   echo "Usage: call with exactly 5 arguments or submit as a batch job with 0 arguments:"
   echo '   ${scr_dir}/compress.csh case_name YYYY-MM-DD-SSSS ensemble_size "sets" "stages"'
   echo '   where '
   echo '   sets   = 1 or more of {clm2 cpl cam cice hist dart} to compress, separated by spaces'
   echo '   stages = 1 or more of stages {input, preassim, postassim, output} to compress.'
   echo ' -OR-'
   echo "   edit compress.csh ; qsub compress.csh"
   exit 17

endif

set cmd = `echo $comp_cmd | cut -d' ' -f1`
if ($cmd == 'gzip') then
   set ext = ''
else if ($cmd == 'gunzip') then
   set ext = '.gz'
else
   echo "ERROR: unrecognized command $cmd.  Don't know which extension to use"
   exit 27
endif

echo "In compress.csh:"
echo "   comp_cmd      = $comp_cmd"
echo "   case_name     = $case_name"
echo "   date          = $ymds"
echo "   ensemble_size = $ensemble_size"
echo "   sets          = $sets"
echo "   stages        = $stages"
echo "   data dir      = $data_dir"

cd $data_dir

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

\rm -f mycmdfile
touch mycmdfile

# 'task' is incremented continuously over all files; components, members, etc.
# 'task' is a running counter of jobs in mycmdfile.
set task = 0

foreach comp ( $sets )
echo "comp = $comp"
switch ($comp)
   # FIXME ... the coupler files may or may not have an instance number in them.
   case {clm2,cpl,cam,cice}:
      set i=1
      while ( $i <= $ensemble_size)
         # E.g. CAM6_80mem.cice_0001.r.2010-07-15-00000.nc
         set file_name = `printf "%s.%s_%04d.r.%s.nc%s" $case_name $comp $i $ymds $ext`
         # echo "   $file_name"
   
         # If the expected file exists, add the compression command 
         if (-f $file_name) then
            @ task++
            echo "$comp_cmd $file_name &> compress_${task}.eo " >> mycmdfile
         # Kluge to get around situations where an earlier job compressed the file,
         # but failed for some other reason, so it's being re-run.
         else if (-f ${file_name}.gz) then
            echo "$file_name already compressed"
         else
            echo 'ERROR: Could not find "'$file_name'" to compress.'
            exit 47
         endif

         @ i++
      end
      breaksw

   case hist:
      # Coupler history (forcing) files, ordered by decreasing size 
      # ha is not a necessary forcing file.  The others can do the job
      # and are much smaller.
      foreach type ( ha2x1d hr2x ha2x3h ha2x1h ha2x1hi )
         # Loop over instance number
         set i=1
         while ( $i <= $ensemble_size)
            # E.g. CAM6_80mem.cpl_0001.ha.2010-07-15-00000.nc
            set file_name = `printf "%s.cpl_%04d.%s.%s.nc%s" $case_name $i $type $ymds $ext`
      
            if (-f $file_name) then
               @ task++
                  echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
            else if (-f ${file_name}.gz) then
               echo "$file_name already compressed"
            else
               echo 'ERROR: Could not find "'$file_name'" to compress.'
               exit 57
            endif

            @ i++
         end
      end
      breaksw

   case dart:
      # It is not worthwhile to compress inflation files ... small, not many files
      # It is also not clear that binary observation sequence files compress effectively.

      # obs_seq.final (no inst)   70% of 1 Gb (ascii) in 35 sec
      # E.g. CAM6_80mem.dart.e.cam_obs_seq_final.2010-07-15-00000
      set file_name = ${case_name}.dart.e.cam_obs_seq_final.${ymds}${ext}
      if (-f $file_name) then
         @ task++
         echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
      endif

      foreach stage ($stages)
         foreach stat ( 'mean' 'sd' )
            # E.g. CAM6_80mem.e.cam_output_mean.2010-07-15-00000.nc
            # E.g. CAM6_80mem.e.cam_output_sd.2010-07-15-00000.nc
            set file_name = ${case_name}.dart.e.cam_${stage}_${stat}.${ymds}.nc${ext}
            if (-f $file_name) then
               @ task++
               echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
            endif
         end

         # Loop over instance number
         set i=1
         while ( $i <= $ensemble_size)
            # E.g. CAM6_80mem.cam_0001.e.preassim.2010-07-15-00000.nc
            set file_name = `printf "%s.cam_%04d.e.%s.%s.nc%s" $case_name $i $stage $ymds $ext`
            if (-f $file_name) then
               @ task++
               echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
            endif
            @ i++
         end
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
   mpiexec_mpt -n $task launch_cf.sh ./mycmdfile

   set mpt_status = $status
   echo "mpt_status = $mpt_status"
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

exit 0


