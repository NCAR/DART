#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# $Id:$

#==========================================================================

# Script to package files found in $DOUT_S_ROOT
# after $reanalysis_cesm/st_archive has sorted them,
# and obs_diag has generated a basic obs space diagnostics NetCDF file.
# The resulting files will be moved to 
#   > a project space for further analysis and use
#   > Campaign Storage for intermediate archiving, 
#     until we want to send them to HPSS (tape).
# Both destinations take time.  They are actually copies.

# >>> Run st_archive and obs_diag before running this script.      <<<
# >>> Check that there's enough disk space.                        <<<
#     It needs 1 Tb more than current usage.                       <<<
#     That's assuming that processing the cpl hist files is first, <<<
#     after which lots of space is freed up.  If it's not,         <<<
#     the lnd history files needs almost 3 Tb additional space.    <<<
# >>> Log in to globus (see mv_to_campaign.csh for instructions).  <<<
# >>> From a casper window (but not 'ssh'ed to data-access.ucar.edu)
#     submit this script from the CESM CASEROOT directory.         <<<

#-----------------------------------------
# Submitting the job to casper (or other NCAR DAV clusters?) requires using slurm.

# Important things to know about slurm:
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submitting a job
# scancel   killing a job
# scontrol  show job <jobID> specifications 
#           (-d for more details, including the script)
#
#==========================================================================
#
#SBATCH --job-name=repack_st_archive
# Output standard output and error to a file named with 
# the job-name and the jobid.
#SBATCH -o %x_%j.eo 
#SBATCH -e %x_%j.eo 
# 80 members
# restarts #SBATCH --ntasks=160 
# forcing files: #SBATCH --ntasks=405 
# partial history set (8 members x 2 types)
#SBATCH --ntasks=16
# 3 members; 
# #SBATCH --ntasks=15 
# #SBATCH --ntasks=1 
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
# #SBATCH --account=P86850054
#SBATCH --account=NCIS0006
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------

# > > > WARNING; cheyenne compute nodes do not have access to Campaign Storage. < < < 
#       Run this on casper using slurm
cd $SLURM_SUBMIT_DIR
# In order to rebase the time variable NCO needs these modules
module load nco gnu udunits

# Needed for mpiexec_mpt:  setenv MPI_SHEPHERD true
setenv date 'date --rfc-3339=ns'

echo "Preamble at "`date`

if (! -f CaseStatus) then
   echo "ERROR: this script must be run from the CESM CASEROOT directory"
   exit 1
endif

setenv CASEROOT      `./xmlquery CASEROOT    --value`
set CASE           = `./xmlquery CASE        --value`
set DOUT_S_ROOT    = `./xmlquery DOUT_S_ROOT --value`
set ensemble_size  = `./xmlquery NINST_ATM   --value`
set line           = `grep '^[ ]*stages_to_write' input.nml`
set stages_all     = (`echo $line[3-$#line] | sed -e "s#[',]# #g"`)

# Non-esp history output which might need to be processed.
# "components" = generic pieces of CESM (used in the archive directory names).
# "models" = component instance names (models, used in file names).
set components     = (lnd  atm ice  rof)
set models         = (clm2 cam cice mosart)

set line = `grep -m 1 save_rest_freq ./assimilate.csh`
set save_rest_freq = $line[4]

if (! -d $DOUT_S_ROOT) then
   echo "ERROR: Missing local archive directory (DOUT_S_ROOT).  "
   echo "       Maybe you need to run st_archive before this script"
   exit 10
endif

# Default mode; archive 5 kinds of data.
# These can be turned off by editing or argument(s).
# do_forcing can only be turned off if archive/cpl/hist/cmds_template exists,
#    or do_history is also turned off.
# Number of tasks required by each section (set according to the max of the 'true's)
# do_forcing     => nens * 5
# do_restarts    => nens + 1
# do_obs_space   => 1
# do_history     => nens * MAX(# history file types.  Currently 2 (CLM))
# do_state_space => 1  (Could be upgraded to use #rest_dates(4-5) * #stats(4))

set do_forcing     = 'true'
# > > > WARNING; if restarts fails when $mm-01 is a Monday, turn off the pre_clean flag,
#                in order to preserve what's in rest/YYYY-MM.
set do_restarts    = 'true'
set do_obs_space   = 'false'
set do_history     = 'true'
set do_state_space = 'true'

# Check whether there is enough project disk space to run this.
# The numbers added to pr_used were harvested from running repack_hwm.csh.
set line = `gladequota | grep ncis0006 `
set pr_used = `echo $line[2] | cut -d'.' -f1`
# Round it up to be safe.
@ pr_used++
if (do_forcing == 'true') then
   @ pr_need = $pr_used + 1
   if ($pr_need > 20) then
      echo "ERROR; not enough project space to run this"
      exit 2
   endif
else if ($do_history == 'true') then
   @ pr_need = $pr_used + 3
   if ($pr_need > 20) then
      echo "ERROR; not enough project space to run this"
      exit 3
   endif
endif

#--------------------------------------------
if ($#argv == 0) then
   # User submitted, independent batch job (not run by another batch job).
   # CASE could be replaced by setup_*, as is done for DART_config.
   # "project" will be created as needed (assuming user has permission to write there).
   # set project    = /glade/p/cisl/dares/Reanalyses/CAM6_2017
   set project    = /glade/p/nsc/ncis0006/Reanalyses
   set campaign   = /gpfs/csfs1/cisl/dares/Reanalyses
   set year  = 2016
# > > > WARNING: if the first day of the month is a Monday,
#       I need to add *_0001.log* files from $archive/logs to rest/YYYY-MM-01-00000
#       and remove the rpointer and .h0. files.
   set month = 12
   set yr_mo = `printf %4d-%02d ${year} ${month}`
   set obs_space  = Diags_NTrS_${yr_mo}

   env | sort | grep SLURM

else if ($#argv == 1) then
   # Request for help; any argument will do.
   echo "Usage:  "
   echo "Before running this script"
   echo "    Run st_archive and obs_diags.csh. "
   echo "    Log in to globus (see mv_to_campaign.csh for instructions)."
   echo "    From a casper window (but not 'ssh'ed to data-access.ucar.edu)"
   echo "    submit this script from the CESM CASEROOT directory. "
   echo "Call by user or script:"
   echo "   repack_st_arch.csh project_dir campaign_dir yr_mo [do_this=false] ... "
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

cd $DOUT_S_ROOT
pwd

# Check that this script has not already run completely for this date.
if ($do_state_space == true && \
    -f ${campaign}/${CASE}/logs/${yr_mo}/da.log.${yr_mo}.tar) then
   echo "ERROR; ${campaign}/${CASE}/logs/${yr_mo}/da.log.${yr_mo}.tar already exists."
   exit 15
endif

if ($do_obs_space != 'true') then
   echo "SKIPPING archiving of obs_space diagnostics"
else if (! -d esp/hist/$obs_space && ! -f esp/hist/${obs_space}.tgz) then
   echo "ERROR: esp/hist/$obs_space does not exist."
   echo "       run obs_diags.csh before this script."
   exit 20
endif

if ($do_forcing == 'false') then
   if (! ( -f cpl/hist/cmds_template || $do_history == 'false')) then
      echo "do_forcing can only be turned off if "
      echo "   archive/cpl/hist/cmds_template exists, "
      echo "   or do_history is also turned off."
      exit 22
   endif
endif


# >>> Remove or update this
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
echo "------------------------"
if ($do_forcing == true) then
   echo "Forcing starts at "`date`
   cd ${DOUT_S_ROOT}/cpl/hist

   # Make a list of the dates (buried in file names).
   set files_dates = `ls ${CASE}.cpl_0001.ha2x1d.${yr_mo}-*.nc*`
   if ($#files_dates == 0) then
      echo "ERROR: there are no ${CASE}.cpl_0001.ha2x1d files.  Set do_forcing = false?"
      exit 23
   endif
   # Separate decompress.csh for each date which needs it.
   foreach d ($files_dates)
      if ($d:e == 'gz') then
         set ymds = $d:r:r:e
         ${CASEROOT}/compress.csh gunzip $CASE $ymds $ensemble_size "hist" "not_used"
         if ($status != 0) then
            echo "ERROR: Compression of coupler history files failed at `date`"
            exit 25
         endif
      endif
   end

   # Make a cmd file to append this month's time series to the yearly file in $project
   # Start with a template of all the instances of one file type.
   if (-f cmds_template) mv cmds_template cmds_template_prev
   touch cmds_template
   set i = 1
   while ($i <= $ensemble_size)
      set NINST = `printf %04d $i`
      set inst_dir = ${project}/${CASE}/cpl/hist/${NINST}
      # "TYPE" will be replaced by `sed` commands below.
      set yearly_file = ${CASE}.cpl_${NINST}.TYPE.${year}.nc

      if (-d $inst_dir) then
         cd ${inst_dir}

         ls  ${CASE}.cpl_${NINST}.*.${year}.nc >& /dev/null
         if ($status == 0) then
            mv   ${CASE}.cpl_${NINST}.*.${year}.nc Previous || exit 28
            # "$init" is a place holder, in the template command, for the existence
            # of a yearly file
            set init = ${inst_dir}/Previous/$yearly_file
         else
            set init = ''
         endif

         cd ${DOUT_S_ROOT}/cpl/hist
     
      else
         mkdir -p $inst_dir
         set init = ''
      endif

      # First try:
      # echo "ncrcat -A -o $yearly_file " \
      #    "${CASE}.cpl_${NINST}.TYPE.${yr_mo}-*.nc &> TYPE_${NINST}.eo " \
      # Apparently casper's ncrcat is not described well by NCO 4.8.1 documentation.
      # 1) It treats -A as (over)writing variables into a file,
      #    rather than appending records to a file.  The latter is done by --rec_apn.
      #    My mistake; the docs do say -A is different from concatenating.
      # 2) If the -o option is used, it does not preserve the time:units attribute 
      #    of the first input file, even though it should.  So the output file 
      #    must be listed after the all of the input files.
      # 3) The output file cannot also be the first input file,
      #    or the time monotonicity may be violated.
      #    This defeats the intent of the "append" mode, but testing confirmed it.

      echo "ncrcat $init  ${DOUT_S_ROOT}/cpl/hist/${CASE}.cpl_${NINST}.TYPE.${yr_mo}*.nc " \
           " ${inst_dir}/$yearly_file &> TYPE_${NINST}.eo " \
           >> cmds_template
      @ i++
   end

   # Append a copy of the template file, modified for each file type, into the command file.
   if (-f mycmdfile) mv mycmdfile mycmdfile_prev
   touch mycmdfile
   @ task = 0
   foreach type (ha2x3h ha2x1h ha2x1hi ha2x1d hr2x)
      sed -e "s#TYPE#$type#g" cmds_template >> mycmdfile
      @ task = $task + $ensemble_size
   end

   echo "   forcing mpirun launch_cf.sh starts at "`date`
   mpirun -n $task ${CASEROOT}/launch_cf.sh ./mycmdfile
   set mpi_status = $status
   echo "   forcing mpirun launch_cf.sh ends at "`date`

   ls *.eo >& /dev/null
   if ($status == 0) then
      grep ncrcat *.eo >& /dev/null
      # grep failure = ncrcat success = "not 0"
      set ncrcat_failed = $status
   else
      # No eo files = failure of something besides g(un)zip.
      echo "cmdfile created no log files for forcing files "
      echo "   and mpi_status of ncrcats = $mpi_status"
      set ncrcat_failed = 0
   endif

   if ($mpi_status == 0 && $ncrcat_failed != 0) then
      rm mycmdfile *.eo $inst_dir:h/*/Previous/*
   else
      echo 'ERROR in repackaging forcing (cpl history) files: See h*.eo, cmds_template, mycmdfile'
      echo '      grep ncrcat *.eo  yielded status '$ncrcat_failed
      exit 50
   endif

   cd ${DOUT_S_ROOT}
endif

# 2) Restart sets, weekly on Monday @ 00Z
#    $DOUT_S_ROOT/rest/YYYY-MM-DD-00000
#    Package by member and date, to allow users 
#      to grab only as many members as needed.
#      a la ./package_restart_members.csh
#    Send to campaign storage using ./mv_to_campaign.csh
# This requires space in $scratch to house the new tar files.
# The files must exist there until globus is done copying them to campaign storage.
# If $scratch has filled with assimilation output, then the forcing file
# section of this script will make enough room for this section.

echo "------------------------"
if ($do_restarts == true) then
   echo "Restarts starts at "`date`
   
   cd ${DOUT_S_ROOT}/rest

   # Pre_clean deals with the feature of mv_to_campaign.csh,
   # which copies all the contents of a directory to campaign storage.
   # If the contents of that directory are left over from a previous repackaging,
   # the directory needs to be cleaned out.
   set pre_clean = true
   # During debugging, it may be helpful to *not* clean out mv_to_campaign's directory.


   # Files_to_save keeps track of whether any restart sets have been packaged
   # for moving to campaign storage.
   set files_to_save = false


   foreach rd (`ls -d ${yr_mo}-*`)
      # Purge restart directories which don't fit in the frequency 
      # defined in setup_* by save_every_Mth_day_restart.
#       echo ' '
#       echo Contents of rest at the start of the dates loop
#       ls -dlt *

      # Prevent archiving files (just do directories).
      if (-f $rd) continue
 
      # The directories we want have only numbers and '-'s.
      echo $rd | grep '[a-zA-Z]' 
      if ($status == 0) continue

      # Ignore directories names that don't have '-'.
      # This doesn't apply to the $yr_mo (mv_to_campaign.csh) directory
      # because it didn't exist when the set of $rd was defined for this loop.
      echo $rd | grep '\-'
      if ($status != 0) continue

      echo ' '
      echo Processing $rd

      set rd_date_parts = `echo $rd | sed -e "s#-# #g"`
      set day_o_month = $rd_date_parts[3]
      set sec_o_day   = $rd_date_parts[4]
      echo year = $year

      if ($pre_clean == true) then
         set pre_clean = false

         # Mv_to_campaign.csh is designed to move everything in a directory to campaign storage.
         # Clean up and/or make a new directory for the repackaged files 
         # and send that directory to mv_to_campaign.csh.
   
         if (-d $yr_mo) then
            rm ${yr_mo}/*
            echo Cleaned out contents of previous restarts from $yr_mo
         else
            mkdir $yr_mo
            echo "Made directory `pwd`/$yr_mo "
            echo "   to store restart members until globus archives them."
         endif
 
      endif

      set purge = 'true'

      # Learn whether save_rest_freq is a string or a number.
      # Character strings must be tested outside of the 'if' statement.
      echo $save_rest_freq | grep '[a-z]'
      if ($status == 0) then
         set purge_date = ${yr_mo}-${day_o_month}
         set weekday = `date --date="$purge_date" +%A`
         if ($weekday == $save_rest_freq && \
             $sec_o_day == '00000') set purge = 'false'
   
      # Numbers can be tested inside the 'if' statement.
      else if (`echo $save_rest_freq | grep '[0-9]'`) then
         if ($day_o_month % $save_rest_freq == 0 && \
             $sec_o_day == '00000') set purge = 'false'

      else
         echo "ERROR: save_every_Mth_day_restart = $save_rest_freq from setup_??? is not supported.  "
         echo "   It must be an integer (< 31) or a day of the week (Monday).  Exiting"
         exit 57 
      endif
   
      if ($purge == 'true') then
         echo "Ignoring restart directory $rd because it doesn't match "
         echo "         save_every_Mth_day_restart = $save_rest_freq and/or $sec_o_day != 00000"
         echo "         It will be removed by purge.csh."

      # This prevents the most recent (uncompressed) restart set from being archived
      # But this may be prevented more reliably by the selection in foreach ( $yr_mo).
      else 
         # YYYY-MM-01-00000 directories are archived by st_archive instead of by assimilate.csh,
         # so they don't have CESM log files in them, which will cause the code below to crash.
         ls ${rd}/*.log.*  >& /dev/null
         if ($status != 0 ) then
            echo "ERROR: $DOUT_S_ROOT/rest/${rd} has no CESM log files."
            echo "       Import them from archive/logs, then remove the rpointer and .h0. files."
            echo "       % ls1 to get the time stamp of the .i., grep for it in logs, "
            echo "       % cp ../../logs/{*_0001,cesm}.log. ... "
            echo "       Then rerun this with do_forcing turned off."
            exit 59
         endif
      
         echo "Exporting restart file set ${rd} "
         echo "   to ${campaign}/${CASE}/rest/${yr_mo} "

         set files_to_save = true

         if (-f mycmdfile) then
            mv mycmdfile mycmdfile_prev
            rm *.eo
         endif
         touch mycmdfile
         set i = 1
         while ($i <= $ensemble_size)
            set NINST = `printf %04d $i`
            echo "tar -c -f ${yr_mo}/${CASE}.${NINST}.alltypes.${rd}.tar "                     \
                           "${rd}/${CASE}.*_${NINST}.*.${rd}.* &>  tar_${NINST}_${rd}.eo "  \
                     "&& rm ${rd}/${CASE}.*_${NINST}.*.${rd}.* &>> tar_${NINST}_${rd}.eo" >> mycmdfile
            @ i++
         end
         # Clean up the rest (non-instance files).
         # 'dart.r' will catch the actual restart files (.rh.)
         # and the restart file (.r.), which is used only by st_archive.
         # Assimilate.csh should always put the .rh. files into the restart sets,
         # so there's no 'if' test around this command (although it might be useful
         # to print a meaningful error message).
         # It's necessary to direct the error output from the rm to the eo file
         # so that it can be examined below and won't necessarily cause an error+exit below.
         echo "tar -c -f ${yr_mo}/${CASE}.infl_log.alltypes.${rd}.tar "     \
                        "${rd}/*{dart.r,log}* &>  tar_inf_log_${rd}.eo "  \
                  "&& rm ${rd}/*{dart.r,log}* &>> tar_inf_log_${rd}.eo " >> mycmdfile
      
         @ tasks = $ensemble_size + 1
         echo "Restart mpirun launch_cf.sh starts at "`date`
         mpirun -n $tasks ${CASEROOT}/launch_cf.sh ./mycmdfile
         set mpi_status = $status
         echo "Restart mpirun launch_cf.sh ends at "`date`" with status "$mpi_status
    
         ls *.eo >& /dev/null
         if ($status == 0) then
            grep tar *.eo | grep -v log >& /dev/null
            # grep failure = tar success = "not 0"
            set ncrcat_failed = $status
         else
            # No eo files = failure of something besides tar.
            set ncrcat_failed = 0
            echo "Restart file set, tar mpi_status = $mpi_status"
         endif
      
         if ($mpi_status == 0 && $ncrcat_failed != 0) then
            rm tar*.eo
         else
            echo 'ERROR in repackaging restart files: See tar*.eo, mycmdfile'
            echo '      grep tar *.eo  yielded status '$ncrcat_failed
            ls -l *.eo
            exit 60
         endif

      endif

      # Copy all the restart tars to campaign storage
      # It's OK to do this within the loop because mv_to_campaign.csh
      # uses the --sync-level argument to globus to move only the new(er) files
      # to campaign storage.
      if ($files_to_save == true) then
         # Remove the empty directory to prevent mv_to_campaign.csh from archiving it.
         rmdir $rd
         if ($status != 0) then
            echo "ERROR; $rd is not empty, Cannot remove it"
            exit 62
         endif
         rm mycmdfile
   
         # Echo the archive command to help with globus error recovery
         # and make it easier to do that on cheyenne as well as on casper.
         echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/rest/$yr_mo " \
                                              " ${campaign}/${CASE}/rest"
         ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/rest/$yr_mo \
                                        ${campaign}/${CASE}/rest

         echo "WARNING: mv_to_campaign.csh ONLY SUBMITS THE REQUEST to globus"
         echo "         The directory must be manually removed after globus completes the copy."
      endif
   end

   cd ${DOUT_S_ROOT}
   
endif


# 3) Obs space diagnostics.
#    Generated separately using ./obs_diags.csh

echo "------------------------"
if ($do_obs_space == true) then
   cd ${DOUT_S_ROOT}/esp/hist
   echo " "
   echo "Location for obs space is `pwd`"
   
   # This is single threaded and takes a long time,
   #  so do_obs_space is usually 'false'
   # and this tar is done in the script that also calls obs_diag: 
   # cesm2_1/diags_rean.csh
   tar -z -c -f ${CASE}.cam_obs_seq_final.${yr_mo}.tgz \
         ${CASE}.dart.e.cam_obs_seq_final.${yr_mo}* &

   # Move the obs space diagnostics to $project.
   set obs_proj_dir = ${project}/${CASE}/esp/hist/${yr_mo}
   if (! -d $obs_proj_dir) mkdir -p $obs_proj_dir

   mv ${obs_space}.tgz $obs_proj_dir
   if ($status == 0) then
      rm -rf $obs_space
   else
      echo "$obs_space.tgz could not be moved.  $obs_space not removed"
   endif

   # Let the tar of obs_seq files finish.
   echo "Waiting for tar of obs_seq files at date - `date --rfc-3339=ns`"
   wait
   if (  -f ${CASE}.cam_obs_seq_final.${yr_mo}.tgz && \
       ! -z ${CASE}.cam_obs_seq_final.${yr_mo}.tgz) then
      mv ${CASE}.cam_obs_seq_final.${yr_mo}.tgz     $obs_proj_dir
      if ($status == 0) then
         rm ${CASE}.dart.e.cam_obs_seq_final.${yr_mo}*
         echo "Moved ${CASE}.cam_obs_seq_final.${yr_mo}.tgz "
         echo "   to $obs_proj_dir"
      endif
   else
      echo "${CASE}.cam_obs_seq_final.${yr_mo}.tgz cannot be moved"
      exit 90
   endif

   cd ${DOUT_S_ROOT}
   
endif

#--------------------------------------------

# 4) CESM history files
#    + CLM output (.h1.) (for Lombardozzi) may need monthly averaging.
#    + .h0. does not., but she only needs daily; concatenate the -00000 files.
#    Should all members be saved?  Then use cmdfile.
#    Or just 1?  Much simpler

echo "------------------------"
if ($do_history == true) then
   # If cam files are too big, do them separately and monthly.
   echo "There are $#components components (models)"
   # Near the beginning of the script:
   # set components     = (lnd  atm ice  rof)
   # set models         = (clm2 cam cice mosart)
   set m = 1
   while ($m <= $#components)
      ls $components[$m]/hist/*h0* >& /dev/null
      if ($status != 0) then
         echo "Skipping $components[$m]/hist"
         @ m++
         continue
      endif
      if ($models[$m] == 'cam' ) then
         ls $components[$m]/hist/*h[^0]* >& /dev/null
         if ($status != 0) then
            echo "Skipping $components[$m]/hist"
            @ m++
            continue
         endif
      endif

      cd $components[$m]/hist
      echo " "
      echo "Location for history is `pwd`"

      set i = 1
      @ comp_ens_size = ( $ensemble_size - $i ) + 1
      while ($i <= $ensemble_size)
         set NINST = `printf %04d $i`
         set inst_dir = ${project}/${CASE}/$components[$m]/hist/${NINST}

         if (-d $inst_dir) then
            cd ${inst_dir}


            # The file form is like yearly_file = ${CASE}.cpl_${NINST}.TYPE.${year}.nc
            # in the forcing file section, but for all TYPEs and a different component.
            ls ${CASE}.$models[$m]_${NINST}.*.${year}.nc >& /dev/null
            if ($status == 0) then
               mkdir -p Previous
               mv ${CASE}.$models[$m]_${NINST}.*.${year}.nc Previous 
            else if ($month != 1) then
               # Exit because if $inst_dir exists there should be a yearly file in it.
               echo "There are no ${CASE}.$models[$m]_${NINST}."'*'".${year}.nc files.  Exiting"
               exit 95
            endif

            cd ${DOUT_S_ROOT}/$components[$m]/hist
         else 
            mkdir -p $inst_dir
         endif

         @ i++
      end
   
      # Make a cmd file to append this month's time series to the yearly file in $project
      # Start with a template of all the instances of one file type.
   
      # This is what's in cmds_template, which is re-used here.
      #       set inst_dir = ${project}/${CASE}/cpl/hist/${NINST}
      #       set yearly_file = ${CASE}.cpl_${NINST}.TYPE.${year}.nc
      #       echo "ncrcat --rec_apn    $yearly_file " \
      #            "${CASE}.cpl_${NINST}.TYPE.${yr_mo}-*.nc ${inst_dir}/$yearly_file &> " \
      #            "TYPE_${NINST}.eo " \
      set cmds_template = cmds_template_$models[$m]

      ls ${DOUT_S_ROOT}/cpl/hist/cmds_template
      if ($status != 0) then
         echo "ERROR: ${DOUT_S_ROOT}/cpl/hist/cmds_template is missing; need it for archiving h# files."
         echo "       It should have been created in section 1 of this script."
         exit 105
      endif
   
      set templ_size = `wc -l ${DOUT_S_ROOT}/cpl/hist/cmds_template`
      if ($templ_size[1] != $comp_ens_size) then
         echo "ERROR: Mismatch of sizes of ${DOUT_S_ROOT}/cpl/hist/cmds_template "
         echo "       and this component's members = $comp_ens_size"
         exit 110
      endif

      sed -e "s#cpl_#$models[$m]_#g;s#cpl#$components[$m]#g"  \
          ${DOUT_S_ROOT}/cpl/hist/cmds_template >! $cmds_template

      # Append a copy of the template file, modified for each file type, into the command file.
      set mycmdfile = mycmdfile_$models[$m]
      if (-f $mycmdfile) mv ${mycmdfile} ${mycmdfile}_prev
      touch ${mycmdfile}

      # The number of history files = SUM(ensemble_size * hist_types_this_comp * dates_this_type)
      @ tasks = 0
      @ type = 0
      while ($type < 10)
         # This learns when there are no more h# types to process.
         # All the desired dates for this type will be appended 
         # to the yearly file for EACH member.
         # Mosart writes out monthly h0 files by default.
         # I don't know whether they have actual monthly averages,
         # or just the last time slot.  There are .rh0 files in $rundir.
         # In any case, they're labeled with YYYY-MM, but no -DD-SSSSS.
         # In contrast to the forcing file list of dates, this list does
         # not include the last "-" because we want to find the Mosart files,
         # but there's no danger of finding the yearly file
         # because it exists only in $project, not locally.
         set dates = `ls *0001.h${type}.${yr_mo}*`
         if ($#dates == 0) break
   
         # There are ensemble_size commands in cmds_template.
         # If cam.h0 ends up with more than PHIS, don't do this if test.
         # and fix the h0 purging in the state_space section.
         if ($models[$m] == 'cam' && $type == 0) then
            # If a single/few cam.h0 files need to be saved (for PHIS):
#             sed -e "s#TYPE#h$type#g" ${cmds_template} | grep _0001 >> ${mycmdfile}
#             @ tasks = $tasks + 1
         else
            sed -e "s#TYPE#h$type#g" ${cmds_template} >> ${mycmdfile}
            @ tasks = $tasks + $ensemble_size
         endif
         @ type++
      end
   
      if (-z $mycmdfile) then
         rm ${cmds_template}  ${mycmdfile} 
         echo "Skipping $components[$m]/hist because $mycmdfile has size 0"
      else
         echo "   history mpirun launch_cf.sh starts at "`date`
         mpirun -n $tasks ${CASEROOT}/launch_cf.sh ./${mycmdfile}
         set mpi_status = $status
         echo "   history mpirun launch_cf.sh ends at "`date`
      
         ls *.eo >& /dev/null
         if ($status == 0) then
            grep ncrcat *.eo >& /dev/null
            # grep failure = ncrcat success = "not 0"
            set ncrcat_failed = $status
         else
            echo "cmdfile created no log files for history files "
            echo "   and mpi_status of ncrcats = $mpi_status"
            set ncrcat_failed = 0
         endif
      
         if ($mpi_status == 0 && $ncrcat_failed != 0) then
            rm ${cmds_template}  ${mycmdfile} *.eo $inst_dir:h/*/Previous/*
         else
            echo "ERROR in repackaging history files: See $components[$m]/hist/"\
                 'h*.eo, cmds_template*, mycmdfile*'
            echo '      grep ncrcat *.eo  yielded status '$ncrcat_failed
            exit 130
         endif
      endif

      cd ${DOUT_S_ROOT}

      @ m++
   end
endif

#--------------------------------------------
# 5) DART diagnostic files: state space
#    + esp/hist/.i.cam_output_{mean,sd}
#    + esp/hist/.e.cam_$stages_{mean,sd}
#    + esp/hist/.rh.cam_$stages_{mean,sd}
#      compress?  85 Mb * 120 dates * 6 uncompressed files= 60 Gb -> 52 Gb/mo.
#                 save 100 Gb / year
echo "------------------------"
if ($do_state_space == true) then
   cd ${DOUT_S_ROOT}/esp/hist
   echo " "
   echo "Location for state space is `pwd`"
   
   mkdir $yr_mo
   foreach stage ($stages_all)
      set ext = e
      if ($stage == output) set ext = i
   
      # Ignoring posterior inflation for now.
      foreach stat  (mean sd priorinf_mean priorinf_sd)
         echo $stat | grep inf 
         if ($status == 0) set ext = rh
         echo "stage, stat, ext = $stage $stat $ext"
   
         if (! -f     ${yr_mo}/${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}.tar) then
            tar -c -f ${yr_mo}/${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}.tar  \
                               ${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}-*
            if ($status == 0) then
               rm ${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}-*
            else
               echo "ERROR: tar -c -f ${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}.tar failed"
               exit 140
            endif
         endif
      end
   end

   # Echo the archive command to help with globus error recovery
   # and make it easier to do that on cheyenne as well as on casper.
   echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/esp/hist/${yr_mo} " \
        " ${campaign}/${CASE}/esp/hist"

   ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/esp/hist/${yr_mo}  \
                                  ${campaign}/${CASE}/esp/hist
 
# The ensemble means are archived every 6 hours, because they're useful and small.
# It also may be useful to have some complete ensembles of model states,
# so those are saved less often (weekly, plus some others sneak in).
# This section archives the ensemble.  
# The members also have a different "file type" than the means
# and are archived to atm/hist, instead of esp/hist
   cd ${DOUT_S_ROOT}/atm/hist
   
   set files = `ls ${CASE}.cam_0001.e.$stages_all[1].${yr_mo}*`
   echo "Files from which atm $stages_all[1] allinst dates will be gathered:"
   echo "  $files"
   if ($#files == 0) then
      echo "There are no .e.$stages_all[1].${yr_mo} files in atm/hist.  Continuing."
   else

      mkdir $yr_mo
      set dates = ()
      foreach f ($files)
         # These files may or may not be compressed, so extracting the date
         # part of the file name is a bit tricky.
         set date = $f:r:e
         if ($date == nc) set date = $f:r:r:e
         set dates = ($dates $date)
      end
      echo " "
      echo "Archiving atm/hist/{$dates} to Campaign Storage"
      
# >>> This could be done a bit faster (1/(#dates*#stages))with mycmdfile code.
      foreach date ($dates)
      foreach stage ($stages_all)
      if ($stage != output) then
   
         tar -c -f ${yr_mo}/${CASE}.cam_allinst.e.${stage}.${date}.tar \
                            ${CASE}.cam_[0-9]*.e.${stage}.${date}*
         if ($status == 0) then
            rm ${CASE}.cam_[0-9]*.e.${stage}.${date}*
         else
            echo "ERROR: tar -c -f ${CASE}.cam_allinst.e.$stage.${date}.tar failed" 
            exit 150
         endif
      endif
      end
      end
      
      # Echo the archive command to help with globus error recovery
      # and make it easier to do that on cheyenne as well as on casper.
      echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/atm/hist/$yr_mo " \
           " ${campaign}/${CASE}/atm/hist"

      ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/atm/hist/$yr_mo \
                                     ${campaign}/${CASE}/atm/hist
   
   endif

   # Archive DART log files (and others?)

   cd ${DOUT_S_ROOT}/logs
   
   # Create a list of files to archive.
   # Logs from more components could be added here.
   set list = ()
   foreach f (`ls da.log*`)
      if ($f:e == 'gz') then
         gunzip $f 
         set f = $f:r
      endif
      grep -l "valid time of model is $year $month" $f >& /dev/null
      if ($status == 0) then
         set list = ($list $f)
      else
         echo "   $f not added to list because 'valid time' $year $month  not found in it. "
      endif
   end

   if ($#list == 0) then
      echo "WARNING: No log files found to archive."
      echo "         da.log list has no files in it."
   endif

   echo "Archiving "
   echo $list

   if (! -d $yr_mo) mkdir $yr_mo
   tar -z -c -f ${yr_mo}/da.log.${yr_mo}.tar $list
   if ($status == 0) then
      rm $list
      # Echo the archive command to help with globus error recovery
      # and make it easier to do that on cheyenne as well as on casper.
      echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/logs/${yr_mo} " \
           "${campaign}/${CASE}/logs"
      ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${DOUT_S_ROOT}/logs/${yr_mo} \
         ${campaign}/${CASE}/logs
   else
      echo "Tar of da.logs of $yr_mo failed.  Not archiving them"
   endif

   cd ${DOUT_S_ROOT}

endif 

#--------------------------------------------

exit 0
#==========================================================================

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
