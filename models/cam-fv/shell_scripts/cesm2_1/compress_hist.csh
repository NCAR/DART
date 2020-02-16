#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Compress just the history files from all of the components,
# one type (h0,...) at a time, which limits the PEs needed
# to the size of the ensemble.

#-----------------------------------------
#SBATCH --job-name=compress_hist
#SBATCH -o %x_%j.eo 
#SBATCH -e %x_%j.eo 
# 80 members
# Each type is done as a separate cmdfile
#SBATCH --ntasks=80 
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
# #SBATCH --account=NCIS0006
# #SBATCH --account=P86850054
#SBATCH --account=R9505652
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------
#PBS  -N compress_hist.csh
#PBS  -A P86850054
#PBS  -q regular
# #PBS  -q share
#PBS  -l select=5:ncpus=36:mpiprocs=36
# #PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=00:10:00
#PBS  -o compress_hist.out
#PBS  -j oe 

# ------------------------------------------------------------------------------
   # Edit these and run as a batch job.
   # 'sets' performs better when ordered by decreasing size (clm2 cpl cam cice hist dart)
   # but this can handle only 1 entry in $sets.
   # -k means keep the input files after creation of the output file.
   # set comp_cmd      = 'gzip -k'
   # No -k option on casper.
   set comp_cmd      = 'gzip'
   set case_name     = f.e21.FHIST_BGC.f09_025.CAM6assim.011
   set CASEROOT      = /glade/work/raeder/Exp/$case_name
   # set ymds          = 2010-07-17-64800
   set ymds          = 2012
   set ensemble_size = 80

   # set data_dir      = ${pr}/${case_name}/cpl/hist
   # set sets          = (cpl)
   # set types         = ( ha2x1d hr2x ha2x3h ha2x1h ha2x1hi )

   # set data_dir      = ${pr}/${case_name}/lnd/hist
   # set data_dir      = /gpfs/csfs1/cisl/dares/Reanalyses/2012_2011SST/lnd/hist
   set data_dir      = /glade/p/nsc/ncis0006/Reanalyses/${case_name}/rof/hist
   set sets          = (mosart)
   set types         = ( h0 )

   # set data_dir      = /glade/scratch/${USER}/${case_name}/run
   # set sets          = (hist dart)
   # set sets          = (clm2 cpl cam cice)
   set stages        = (none)

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
         while ( $i <= $ensemble_size)
            # E.g. CAM6_80mem.cice_0001.r.2010-07-15-00000.nc
            set file_name = \
                `printf "%04d/%s.%s_%04d.%s.%s.nc%s" $i $case_name $comp $i ${type} $ymds $ext`
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
         mpiexec_mpt -n $task ${CASEROOT}/launch_cf.sh ./mycmdfile
         set mpi_status = $status
      else if ($?SLURM_SUBMIT_DIR) then
         mpirun      -n $task ${CASEROOT}/launch_cf.sh ./mycmdfile
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
