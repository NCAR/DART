#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# DART $Id$

#==========================================================================

# >>> Run st_archive and obs_diags.csh before running this script. <<<
# >>> Edit compress.csh to make it decompress the (cpl/hist) files,
#     in cases where they are compressed (gzipped).
# >>> Log in to globus
# >>> Run this script from the CESM CASEROOT directory. <<<

# Script to package files found in $DOUT_S_ROOT
# after $reanalysis_cesm/st_archive has sorted them,
# and obs_diags.csh has generated basic obs space diagnostics.
# The resulting files will be moved to 
#   > a project space for further analysis and use
#   > Campaign Storage for temporary archiving, 
#     until we want to send them to HPSS (tape).
# Both destinations take time.  They are actually copies.

#-----------------------------------------
# Submitting the job to casper (or other NCAR DAV clusters?) requires using slurm.

# Important things to know about slurm:
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submitting a job
# scancel   killing a job
#
#==========================================================================
#
#SBATCH --job-name=repack_st_archive
#SBATCH --ntasks=1 
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
#SBATCH --account=P86850054
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------
#PBS  -N repack_st_archive
#PBS  -A P86850054
# #PBS  -q premium
#PBS  -q share
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
# Request enough processors for command file commands to handle the larger of 
#    > ensemble size x 5 (forcing file types), 
#    > (ensemble size + 1(mean_sd_log)) x # dates being saved (4-5 / month.  )
# 81->405 #PBS  -l select=12:ncpus=34:mpiprocs=34
#PBS  -l select=1:ncpus=15:mpiprocs=15
#PBS  -l walltime=02:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M raeder@ucar.edu
# Send standard output and error to this file.
# It's helpful to use the $CASE name here.
#PBS  -o repack_st_archive.eo
#PBS  -j oe 
#--------------------------------------------

if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
else if ($?SLURM_SUBMIT_DIR) then
   cd $SLURM_SUBMIT_DIR
endif

setenv MPI_SHEPHERD true
setenv date 'date --rfc-3339=ns'

echo "Preamble at "`date`

if (! -f CaseStatus) then
   echo "ERROR: this script must be run from the CESM CASEROOT directory"
   exit 1
endif

set CASEROOT       = $cwd
set CASE           = $CASEROOT:t
set local_arch     = `./xmlquery DOUT_S_ROOT --value`
set ensemble_size  = `./xmlquery NINST_ATM   --value`

set line = `grep -m 1 save_every_Mth_day_restart ./setup*`
set save_rest_freq = $line[3]


if (! -d $local_arch) then
   echo "ERROR: Missing local_arch.  "
   echo "       Maybe you need to run st_archive before this script"
   exit 10
endif

# Default mode; archive 4 kinds of data.
# These can be turned off by editing or argument(s).
set do_forcing     = 'true'
set do_restarts    = 'true'
set do_obs_space   = 'false'
set do_state_space = 'false'

#--------------------------------------------
if ($#argv == 0) then
   # User submitted, independent batch job (not run by another batch job).
   # CASE could be replaced by setup_*, as is done for DART_config.
   # "project" will be created as needed (assuming user has permission to write there).
   # set project    = /glade/p/cisl/dares/Reanalyses/CAM6_2017
   set project    = /glade/p/nsc/ncis0006/Reanalyses/CAM6_2017
   set campaign   = /gpfs/csfs1/cisl/dares/Reanalyses/CAM6_2017
   set yr_mo      = '2017-01'
   set obs_space  = Diags_NTrS_${yr_mo}
   # set DART       = ~/DART/reanalysis

   env | sort | grep SLURM

else if ($#argv == 1) then
   # Request for help; any argument will do.
   echo "Usage: call by user or script:"
   echo "   repack_st_archive.csh project_dir campaign_dir yr_mo [do_this=false] ... "
   echo "      project_dir    = directory where $CASE.dart.e.cam_obs_seq_final.$date.nc are"
   echo "      campaign_dir   = directory where $CASE.dart.e.cam_obs_seq_final.$date.nc are"
   echo "      yr_mo =        = Year and month to be archived, in form YYYY-MO."
   echo "      do_this=false  = Turn off one (or more) of the archiving sections."
   echo "                       'this' = {forcing,restarts,obs_space,state_space}."
   echo "                       No quotes, no spaces."
   exit

else
   # Script run by another (batch) script or interactively.
   set project  = $1
   set campaign = $2
   set yr_mo    = $3
   set obs_space  = Diags_NTrS_${yr_mo}
   # These arguments to turn off parts of the archiving must have the form do_obs_space=false ;
   # no quotes, no spaces.
   if ($?4) set $4
   if ($?5) set $5
   if ($?6) set $6
endif

cd $local_arch
pwd

if ($do_obs_space != 'true') then
   echo "SKIPPING archiving of obs_space diagnostics"
else if (! -d esp/hist/$obs_space) then
   echo "ERROR: esp/hist/$obs_space does not exist."
   echo "       run obs_diags.csh before this script."
   exit 20
endif



#==========================================================================
# Where to find files, and what to do with them
#==========================================================================
# 1) Forcing files
#    $DOUT_S_ROOT/cpl/hist
#    Eventually want(?) each instance's filetype to have a year of data in it.
#       (max 8.2 Gb, seems managable)
#       So no packaging at this point?  Or would monthly be useful?
# ?  Lots of small files is not really a problem when moving to $project?
#      5*80*4*365= 600,000/year, then combined down to 5*80 = 400
#      Could ncrcat each month's onto the existing yearly files.
#      (decompress? and) ncrcat into monthly files?
#    Stream files are specified by directory and file name,
# ?    with variables(?) like NINST and RUNYEAR being filled in at run time.
#      /glade/p/cisl/dares/Reanalyses/CAM6_RUNYEAR/archive/cpl/hist/NINST
#      ${rean_name}.cpl_NINST.${filetype}.RUNYEAR.nc     (5 cpl hist filetypes)
#    Package by member, to allow users to grab only as many members as needed.
#      a la ./package_restart_members.csh
#    Move to $project 
#       This will make room for the next, restart file section to operate.
echo "Forcing starts at "`date`
if ($do_forcing == true) then
   cd cpl/hist

   # Make a list of the dates (buried in file names).
   set files_dates = `ls ${CASE}.cpl_0001.ha2x1d.${yr_mo}-*.nc*`
   if ($#files_dates == 0) then
      echo "ERROR: there are no ${CASE}.cpl_0001.ha2x1d files.  Set do_forcing = false?"
      exit
   else if ($files_dates[1]:e == 'gz') then
      # Separate (de)compress.csh for each date
      foreach d ($files_dates)
         set ymds = $d:r:r:e
         ${CASEROOT}/compress.csh $CASE $ymds $ensemble_size "hist" "not_used"
         if ($status != 0) then
            echo "ERROR: Compression of coupler history files and DART files failed at `date`"
            exit 25
         endif
      end
   endif

   # Make a cmd file to append this month's time series to the yearly file in $project
   # Start with a template of all the instances of one file type.
   set year = `echo $yr_mo | cut -d'-' -f1`
   if (-f tmp_cmd) mv tmp_cmd tmp_cmd_prev
   touch tmp_cmd
   set i = 1
   while ($i <= $ensemble_size)
      set NINST = `printf %04d $i`
      set inst_dir = ${project}/${CASE}/cpl/hist/${NINST}
      if (! -d $inst_dir) mkdir -p $inst_dir

      set yearly_file = ${inst_dir}/${CASE}.cpl_${NINST}.TYPE.${year}.nc

      # && is bash for 'do the next command if the previous succeeded.
      echo "ncrcat -A -o $yearly_file   ${CASE}.cpl_${NINST}.TYPE.*.nc &> TYPE_${NINST}.eo " \
                                 "&& rm ${CASE}.cpl_${NINST}.TYPE.*.nc " \
           >> tmp_cmd

      @ i++
   end

   # Append a copy of the template file, modified for each file type, into the command file.
   if (-f mycmdfile) mv mycmdfile mycmdfile_prev
   touch mycmdfile
   @ task = 0
   foreach type (ha2x3h ha2x1h ha2x1hi ha2x1d hr2x)
      sed -e "s#TYPE#$type#g" tmp_cmd >> mycmdfile
      @ task = $task + $ensemble_size
   end

   echo "   mpiexec_mpt starts at "`date`
   mpiexec_mpt -n $task launch_cf.sh ./mycmdfile
   set mpt_status = $status
   echo "   mpiexec_mpt ends at "`date`

   ls *.eo
   if ($status == 0) then
      grep ncrcat *.eo >& /dev/null
      # grep failure = ncrcat success = "not 0"
      set gr_stat = $status
   else
      # No eo files = failure of something besides g(un)zip.
      set gr_stat = 0
      echo "Forcing file ncrcat mpt_status = $mpt_status"
   endif

   if ($mpt_status == 0 && $gr_stat != 0) then
      rm tmp_cmd mycmdfile *.eo
   else
      echo 'ERROR in repackaging forcing (cpl history) files: See h*.eo, tmp_cmd, mycmdfile'
      echo '      grep ncrcat *.eo  yielded status '$gr_stat
      exit 20
   endif

   cd ../..
endif

# 2) Restart sets, weekly on Monday @ 00Z
#    $DOUT_S_ROOT/rest/YYYY-MM-DD-00000
#    Package by member and date, to allow users 
#      to grab only as many members as needed.
#      a la ./package_restart_members.csh
#    Send to campaign storage using ./mv_to_campaign.csh
#    >>> Needs 
# This requires space in $scratch to house the new tar files.
# The files must exist there until globus is done copying them to campaign storage.
# If $scratch has filled with assimilation output, then the first, forcing file,
# section of this script will make enough room for this section.

echo "Restarts starts at "`date`
if ($do_restarts == true) then
   
   cd rest

   set pre_clean = true
   foreach rd (`ls -d *`)
      # Purge restart directories which don't fit in the frequency 
      # defined in setup_* by save_every_Mth_day_restart.
      echo ' '
      echo Contents of rest at the start of the dates loop
      ls -dlt *

      echo ' '
      echo Processing $rd
      if (-f $rd) continue

      set rd_date_parts = `echo $rd | sed -e "s#-# #g"`
      set year        = $rd_date_parts[1]
      echo year = $year
   
      echo $rd | grep '[a-zA-Z]' 
      if ($status == 0) continue

      echo $rd | grep '\-'
      if ($status == 0) then
         set month       = $rd_date_parts[2]
         set day_o_month = $rd_date_parts[3]
         set sec_o_day   = $rd_date_parts[4]
      else
         echo continuing
         continue
      endif

      if ($pre_clean == true) then
         set pre_clean = false

         # Mv_to_campaign.csh is designed to move everything in a directory to campaign storage.
         # Clean up and/or make a new directory for the repackaged files 
         # and send that directory to mv_to_campaign.csh.
         @ prev_year = $year - 1
         echo prev_year = $prev_year
         if (-d $prev_year) rm -rf $prev_year
   
         if (-d $year) then
            rm $year/*
            echo Cleaned out contents of previous restarts from $year
         else
            mkdir $year
            echo Made directory `pwd`/$year to store restart members until globus archives them.
            # Date for use in mv_to_campaign.csh
            set yr_mo = ${year}-${month}
         endif
   
      endif
      # Learn whether save_rest_freq is a string or a number.

      # Character strings must be tested outside of the 'if' statement.
      set purge = 'true'

      echo $save_rest_freq | grep '[a-z]'
      if ($status == 0) then
         set purge_date = $rd_date_parts[1]-$rd_date_parts[2]-$rd_date_parts[3]
         set weekday = `date --date="$purge_date" +%A`
         if ($weekday == $save_rest_freq) set purge = 'false'
   
      # Numbers can be tested inside the 'if' statement.
      else if (`echo $save_rest_freq | grep '[0-9]'`) then
         if ($day_o_month % $save_rest_freq == 0 && $sec_o_day == '00000') set purge = 'false'

      endif
   
      # Don't purge or archive restart sets at the beginning of the month,
      # which may be needed for the next month.
      set afile = `ls $rd/*cam_0001.r.*${month}-${day_o_month}*`
   
      if ($purge == 'true') then
         echo "Removing restart directory $rd because it doesn't match "
         echo "         save_every_Mth_day_restart = $save_rest_freq and/or $sec_o_day != 00000"
         rm -rf $rd

      # This prevents the most recent (uncompressed) restart set from being archived
      else if ($afile:e == 'gz') then
         echo "Exporting restart file set ${rd} to ${campaign}/${CASE}/rest/${year} "
         echo "WARNING: this command ONLY SUBMITS THE REQUEST to globus"
         echo "         The directory must be manually removed after globus completes the copy."

         if (-f mycmdfile) then
            mv mycmdfile mycmdfile_prev
            rm *.eo
         endif
         touch mycmdfile
         set i = 1
         while ($i <= $ensemble_size)
            set NINST = `printf %04d $i`
            echo "tar -c -f $year/${CASE}.${NINST}.alltypes.${rd}.tar "                     \
                           "${rd}/${CASE}.*_${NINST}.*.${rd}.* &> tar_${NINST}_${rd}.eo "        \
                     "&& rm ${rd}/${CASE}.*_${NINST}.*.${rd}.* " >> mycmdfile
            @ i++
         end
         # Clean up the rest (non-instance files).
         # 'dart.r' will catch the actual restart files (.rh.)
         # and the restart file (.r.), which is used only by st_archive.
         # Assimilate.csh should always put the .rh. files into the restart sets,
         # so there's no 'if' test around this command (although it might be useful
         # to print a meaningful error message).
         echo "tar -c -f $year/${CASE}.infl_log.alltypes.${rd}.tar "                     \
                        "${rd}/*{dart.r,log}* &> tar_inf_log_${rd}.eo "        \
                  "&& rm ${rd}/*{log}* " >> mycmdfile
      
         echo "   mpiexec_mpt starts at "`date`
         @ tasks = $ensemble_size + 1
         mpiexec_mpt -n $tasks launch_cf.sh ./mycmdfile
         set mpt_status = $status
         echo "   mpiexec_mpt ends at "`date`
      
         ls *.eo
         if ($status == 0) then
            grep tar *.eo >& /dev/null
            # grep failure = tar success = "not 0"
            set gr_stat = $status
         else
            # No eo files = failure of something besides tar.
            set gr_stat = 0
            echo "Restart file set, tar mpt_status = $mpt_status"
         endif
      
         if ($mpt_status == 0 && $gr_stat != 0) then
            rm tar*.eo
         else
            echo 'ERROR in repackaging restart files: See tar*.eo, mycmdfile'
            echo '      grep tar *.eo  yielded status '$gr_stat
            exit 20
         endif
      
      else
         echo Did not remove or tar $rd
         
      endif
   end
   # Copy all the restart tars to campaign storage
   # 2019-4-26; Can't ever do from cheyenne batch nodes.
   #            Casper is down this week.
   #            Print the command so that I can run it interactively after this job is done.
   #                          case time_str local_dir                destination
   # ${CASEROOT}/mv_to_campaign.csh $CASE ${year}-${month} ${local_arch}/rest/$year ${campaign}/${CASE}/rest/$year
   echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/rest/$year " \
        " ${campaign}/${CASE}/rest/$year"
   cd ..
endif


# 3) Obs space diagnostics.
#    Generated separately using ./obs_diags.csh
#    Move to $project

# 4) DART diagnostic files
#    + esp/hist/.i.cam_output_{mean,sd}
#    + esp/hist/.e.cam_$stages_{mean,sd}
#    + atm/hist/$inst.e.$stages
#      compress?  85 Mb * 80 * #_dates = 816 Gb -> 710 Gb/mo.
#      Maybe not worth it.  The savings is over a Tb/year.
#    Send to campaign storage using ./mv_to_campaign.csh

# >>> Run this job;
#     + as a batch job, as well as interactive.
#       - PBS
#       - slurm?
#       - args to specify which things to archive, casename(?), pathnames(?)
#     + on casper?  To take advantage of matlab.  
#       NO need; obs_diags.csh run before this.
#     + with switches to turn on individual parts 
#       (to more easily recover from failures).
#     + with echoes about progress and failures.
#     + using cmdfile to do parts simultaneously
#        [Save this worry for the future package_reanalysis.csh
#         - obs_diag + matlab will take the longest; 1 task]
#        - compress obs_seq.finals (and other 'dart' files in ./compress.csh)
#          after obs_diag, parallel to matlab.
#        - tar month of related files into useful units
#            > (compressed) obs_seq.finals
#            > .e.cam_$stages_{mean,sd}
#            > .i.cam_output_{mean,sd}
#            > log files?
#            > .dart.rh.*inf* files (all stages? or just output?)
#        - compressing and/or tarring restarts can take a while; tasks 2,...?
#        - 
#     + after each month of assimilation is done.

exit
#==========================================================================
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
                                                              
