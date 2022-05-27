#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#==========================================================================

# Script to package yearly files found in $data_proj_space 
# (e.g. /glade/p/nsc/ncis0006/Reanalyses/f.e21.FHIST_BGC.f09_025.CAM6assim.011)
# after repack_st_arch.csh has created them,
# and the matlab scripts have generated the obs space pictures. 
# The resulting files will be moved to 
#   > Campaign Storage for intermediate archiving, 
#     until we want to send them to the RDA.
# This takes time; it actually copies.

# FIXME KEVIN - describe how this is supposed to be run/used.
# >>> Run repack_st_arch.csh before running this script. <<<
# >>> Log in to globus (see mv_to_campaign.csh for instructions).
# >>> Purge extraneous files from $data_project_space/$data_CASE/esp/hist, 
#     since (the new parts of) that whole directory will be archived.
# >>> From a casper window (but not 'ssh'ed to data-access.ucar.edu)
#     submit this script from the CESM CASEROOT directory. <<<


#==========================================================================

# Important things to know about slurm:
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submitting a job
# scancel   killing a job
# scontrol  show job <jobID> specifications 
#           (-d for more details, including the script)

#==========================================================================

#SBATCH --job-name=repack_project
# Output standard output and error to a file named with 
# the job-name and the jobid.
#SBATCH -o %x_%j.eo 
#SBATCH -e %x_%j.eo 
# 80 members (1 type at a time)
#SBATCH --ntasks=80 
#SBATCH --time=04:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------
#PBS  -N repack_project.csh
#PBS  -A P86850054
#PBS  -q casper
#PBS  -l select=3:ncpus=36:mpiprocs=36
#PBS  -l walltime=04:00:00
# Make output consistent with SLURM (2011-2019 files)
#PBS  -o repack_project_2020+COSMIC2lnd.eo
#PBS  -j oe 
#PBS  -k eod
#-----------------------------------------

if ($?SLURM_SUBMIT_DIR) then
   cd $SLURM_SUBMIT_DIR
   env | sort | grep SLURM
else if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
   env | sort | grep PBS
endif

# Needed for globus+ncar_py (but only in mv_to_campaign.csh?), but not for gci.
# module load nco gnu 

# Needed for mpiexec_mpt:  
# setenv MPI_SHEPHERD true
setenv date_rfc 'date --rfc-3339=ns'

echo "Preamble at "`$date_rfc`

if (! -f CaseStatus) then
   echo "ERROR: this script must be run from the CESM CASEROOT directory"
   exit 1
endif

# Get CASE environment variables from the central variables file.
source ./data_scripts.csh
echo "data_CASE     = $data_CASE"
echo "data_NINST    = $data_NINST"
echo "data_year     = $data_year"
echo "data_CASEROOT   = $data_CASEROOT"
echo "data_proj_space = $data_proj_space"
echo "data_campaign   = ${data_campaign}"

# Non-esp history output which might need to be processed.
# "components" = generic pieces of CESM (used in the archive directory names).
# "models" = component instance names (models, used in file names).
# "cpl" needs to be first (as in repack_st_arch.csh) so that the cmdfile template
# is created there and can be found by other components.
set components     = (cpl lnd  atm ice  rof)
set models         = (cpl clm2 cam cice mosart)
# set components     = (rof)
# set models         = (mosart)

# Default mode; archive 2 kinds of data.
# These can be turned off by editing or argument(s).
# Number of tasks required by each section (request in the slurm directives 
#     according to the max of the 'true's)
# set do_obs_space   = 1       (This doesn't take long; the whole esp/hist directory is sent to mv_to_campaign.csh)
# set do_history     = nens    (each file type is a separate, 80 task, cmdfile command)
# do_forcing is not needed because it can be handled the same as the 
# {other components}/hist directories.

set do_obs_space   = 'true'
set do_history     = 'true'

#--------------------------------------------
if ($#argv != 0) then
   # Request for help; any argument will do.
   echo "Usage:  "
   echo "Before running this script"
   echo "    Run repack_st_archive.csh. "
   echo "Batch job"
   echo "    submit this script from the CESM CASEROOT directory. "
#    echo "Call by user or script:"
#    echo "   repack_project.csh project_dir campaign_dir [do_this=false] ... "
#    echo "      project_dir    = directory where the yearly files are accumulated by repack_st_arch.csh"
#    echo "      campaign_dir   = directory where compressed yearly files and obs space files will be copied"
#    echo "      do_this=false  = Turn off one (or more) of the archiving sections."
#    echo "                       'this' = {obs_space,hist}."
#    echo "                       No quotes, no spaces."
   exit
endif

# User submitted, independent batch job (not run by another batch job).
# CASE could be replaced by setup_*, as is done for DART_config.
# "data_proj_space" will be created as needed (assuming user has permission to write there).

#==========================================================================
# Where to find files, and what to do with them
#==========================================================================
# 1) Obs space diagnostics.
#    Generated separately using diags_rean.csh and matlab scripts

echo "------------------------"
if ($do_obs_space == true) then
   cd ${data_proj_space}/${data_CASE}/esp/hist
   echo " "
   echo "Location for obs space is `pwd`"

   # Obs space files are already compressed.

   # mv_to_campaign.csh uses 'transfer --sync_level mtime'
   # so the whole esp/hist directory can be specified,
   # but only the files which are newer than the CS versions will be transferred.

   ${data_CASEROOT}/mv_to_campaign.csh \
      $data_year ${data_proj_space}/${data_CASE}/esp/hist/ \
      ${data_campaign}/${data_CASE}/esp/hist

   # If this is successful, it is safe to remove the original obs_seq_final files
   # from $DOUT_S_ROOT/esp/hist, since there will be tarred versions on $project
   # and CS.
   # Do that here?  Or manually?
   
endif

cd ${data_proj_space}/${data_CASE}

#--------------------------------------------

# 2) CESM history files

echo "------------------------"
if ($do_history == true) then
   echo "$#components components (models) were requested"
   echo " "
   set m = 1
   while ($m <= $#components)
      if (! -d $components[$m]/hist) then
         echo "Skipping $components[$m] because there are no history files"
         @ m++
         continue
      endif 

      cd ${data_proj_space}/${data_CASE}/$components[$m]/hist

      echo "============================"
      echo "Location for history is `pwd`"

      if ($components[$m] == 'cpl') then
         set types = ( ha2x1d hr2x ha2x3h ha2x1h ha2x1hi )
      else
         ls 0001/*h0* >& /dev/null
         if ($status != 0) then
            echo "Skipping $components[$m]/hist"
            @ m++
            continue
         endif

         set types = ()
         set n = 0
         while ($n < 10)
            ls 0001/*h${n}* 
            if ($status != 0) then
               @ n = $n - 1
               break
            endif

            set types = ($types h$n)
            @ n++
         end
      endif

      set t = 1
      while ($t <= $#types)
         if ($models[$m] == 'cam' && $t == 1) then
            # If cam.h0 ends up with more than PHIS, comment this out
            # and fix the h0 purging in the state_space section.
            # Actually; skip cam*.h0. because of the purging done by assimilate.csh.
#             sed -e "s#TYPE#h$type#g" ${cmds_template} | grep _0001 >> ${mycmdfile}
#             @ tasks = $tasks + 1
            @ t++
            continue
         else
            echo "   ----------------------"
            echo "   Processing $models[$m] $types[$t]"
         endif

         # Make a cmd file to compress this year's history file(s) in $data_proj_space.
         if (-f cmdfile) mv cmdfile cmdfile_prev
         touch cmdfile

         set tasks = 0
         set i = 1
         while ($i <= $data_NINST)
            set inst = `printf %04d $i`
            set yearly_file = ${data_CASE}.$models[$m]_${inst}.$types[$t].${data_year}.nc

            if (-f ${inst}/${yearly_file}) then
               echo "gzip ${inst}/${yearly_file} &> $types[$t]_${inst}.eo " >> cmdfile
               @ tasks++
            endif
            @ i++
         end

         if (-z cmdfile) then
            echo "WARNING: cmdfile has size 0, hopefully because type $types[$t] was already done"
            @ t++
            continue
         endif

         echo "   history mpirun launch_cf.sh of compression starts at "`$date_rfc`
         mpirun -n $tasks ${data_CASEROOT}/launch_cf.sh ./cmdfile
         set mpi_status = $status
         echo "   history mpirun launch_cf.sh ends at "`$date_rfc`
      
         ls *.eo > /dev/null
         if ($status == 0) then
            grep gzip *.eo >& /dev/null
            # grep failure = gzip success = "not 0"
            set gr_stat = $status
         else
            echo "No eo files = failure of something besides g(un)zip."
            echo "   History file gzip mpi_status = $mpi_status"
            set gr_stat = 0
         endif
      
         if ($mpi_status == 0 && $gr_stat != 0) then
            rm cmdfile *.eo
         else
            echo "ERROR in repackaging history files: See $components[$m]/hist/"\
                 'h*.eo, cmdfile'
            echo '      grep gzip *.eo  yielded status '$gr_stat
            exit 130
         endif

         @ t++
      end

      echo "   Calling mv_to_campaign.csh to copy compressed files"
      ${data_CASEROOT}/mv_to_campaign.csh \
         $data_year ${data_proj_space}/${data_CASE}/$components[$m]/hist/  \
         ${data_campaign}/${data_CASE}/$components[$m]/hist
      echo "   Done with mv_to_campaign.csh"
      echo " "
 
      cd ${data_proj_space}

      @ m++
   end
endif

exit 0
