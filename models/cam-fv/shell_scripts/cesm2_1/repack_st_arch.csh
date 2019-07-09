#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# DART $Id$

#==========================================================================

# Script to package files found in $DOUT_S_ROOT
# after $reanalysis_cesm/st_archive has sorted them,
# and obs_diags.csh has generated basic obs space diagnostics.
# The resulting files will be moved to 
#   > a project space for further analysis and use
#   > Campaign Storage for temporary archiving, 
#     until we want to send them to HPSS (tape).
# Both destinations take time.  They are actually copies.

# >>> Run st_archive and obs_diags.csh before running this script. <<<
# >>> Log in to globus (see mv_to_campaign.csh for instructions).
# >>> From a casper window (but not 'ssh'ed to data-access.ucar.edu)
#     submit this script from the CESM CASEROOT directory. <<<

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
# 80 members
# #SBATCH --ntasks=405 
# 3 members; 
#SBATCH --ntasks=15 
# #SBATCH --ntasks=1 
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
#SBATCH --account=NCIS0006
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------

# > > > WARNING; cheyenne compute nodes do not have access to Campaign Storage. < < < 
#       Run this on casper using slurm


# #PBS  -N repack_st_archive
# #PBS  -A NCIS0006
# # Resources I want:
# #    select=#nodes
# #    ncpus=#CPUs/node
# #    mpiprocs=#MPI_tasks/node
# # Request enough processors for command file commands to handle the larger of 
# #    > ensemble size x 5 (forcing file types), 
# #    > (ensemble size + 1(mean_sd_log)) x # dates being saved (4-5 / month.  )
# #    > ensemble size * max(CAM h#, CLM h#, ...)
# # 81->405 
# 
# #PBS  -q regular
# # #PBS  -l select=12:ncpus=34:mpiprocs=34
# # 3 members:
# #PBS  -l select=1:ncpus=15:mpiprocs=15
# 
# # #PBS  -q share
# # #PBS  -l select=1:ncpus=1:mpiprocs=1
# 
# #PBS  -l walltime=02:00:00
# # Send email after a(bort) or e(nd)
# #PBS  -m ae
# #PBS  -M raeder@ucar.edu
# # Send standard output and error to this file.
# # It's helpful to use the $CASE name here.
# #PBS  -o Test2.eo
# #PBS  -j oe 
# #--------------------------------------------

# if ($?PBS_O_WORKDIR) then
#    cd $PBS_O_WORKDIR
# else if ($?SLURM_SUBMIT_DIR) then
cd $SLURM_SUBMIT_DIR
# In order to rebase the time variable NCO needs these modules
module load nco gnu udunits
# endif

# Needed for mpiexec_mpt:  setenv MPI_SHEPHERD true
setenv date 'date --rfc-3339=ns'

echo "Preamble at "`date`

if (! -f CaseStatus) then
   echo "ERROR: this script must be run from the CESM CASEROOT directory"
   exit 1
endif

setenv CASEROOT $cwd
set CASE           = $CASEROOT:t
set local_arch     = `./xmlquery DOUT_S_ROOT --value`
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

if (! -d $local_arch) then
   echo "ERROR: Missing local_arch.  "
   echo "       Maybe you need to run st_archive before this script"
   exit 10
endif

# Default mode; archive 5 kinds of data.
# These can be turned off by editing or argument(s).
# do_forcing can only be turned off if archive/cpl/hist/cmds_template exists,
#    or do_history is also turned off.
# Number of tasks required by each section (set according to the max of the 'true's)
# set do_forcing     = nens * 5
# set do_restarts    = nens + 1
# set do_obs_space   = 1
# set do_history     = nens * MAX(# history file types.  Currently 2 (CLM))
# set do_state_space = 1  (Could be upgraded to use #rest_dates(4-5) * #stats(4))

set do_forcing     = 'true'
set do_restarts    = 'true'
set do_obs_space   = 'true'
set do_history     = 'true'
set do_state_space = 'true'

#--------------------------------------------
if ($#argv == 0) then
   # User submitted, independent batch job (not run by another batch job).
   # CASE could be replaced by setup_*, as is done for DART_config.
   # "project" will be created as needed (assuming user has permission to write there).
   # set project    = /glade/p/cisl/dares/Reanalyses/CAM6_2017
   set project    = /glade/p/nsc/ncis0006/Reanalyses/CAM6_2010
   set campaign   = /gpfs/csfs1/cisl/dares/Reanalyses/CAM6_2010
   set year  = 2011
   set month = 1
   set yr_mo = `printf %4d-%02d ${year} ${month}`
   set obs_space  = Diags_NTrS_${yr_mo}
   # set DART       = ~/DART/reanalysis

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
   cd cpl/hist

   # Make a list of the dates (buried in file names).
   set files_dates = `ls ${CASE}.cpl_0001.ha2x1d.${yr_mo}-*.nc*`
   if ($#files_dates == 0) then
      echo "ERROR: there are no ${CASE}.cpl_0001.ha2x1d files.  Set do_forcing = false?"
      exit
   else if ($files_dates[1]:e == 'gz') then
      # Separate decompress.csh for each date
      foreach d ($files_dates)
         set ymds = $d:r:r:e
         ${CASEROOT}/compress.csh gunzip $CASE $ymds $ensemble_size "hist" "not_used"
         if ($status != 0) then
            echo "ERROR: Compression of coupler history files and DART files failed at `date`"
            exit 25
         endif
      end
   endif

   # Make a cmd file to append this month's time series to the yearly file in $project
   # Start with a template of all the instances of one file type.
   # Not needed since year is defined at the start
   # set year = `echo $yr_mo | cut -d'-' -f1`
   if (-f cmds_template) mv cmds_template cmds_template_prev
   touch cmds_template
   set i = 1
   while ($i <= $ensemble_size)
      set NINST = `printf %04d $i`
      set inst_dir = ${project}/${CASE}/cpl/hist/${NINST}
      if (! -d $inst_dir) mkdir -p $inst_dir

      set yearly_file = ${CASE}.cpl_${NINST}.TYPE.${year}.nc

      # && is bash for 'do the next command if the previous succeeded.
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
      #
      if (-f ${inst_dir}/$yearly_file) then
         mv ${inst_dir}/$yearly_file . 
         set init = $yearly_file
      else
         set init = ''
      endif
      echo "ncrcat --rec_apn  $init ${CASE}.cpl_${NINST}.TYPE.${yr_mo}-*.nc " \
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

   ls *.eo > /dev/null
   if ($status == 0) then
      grep ncrcat *.eo >& /dev/null
      # grep failure = ncrcat success = "not 0"
      set gr_stat = $status
   else
      # No eo files = failure of something besides g(un)zip.
      set gr_stat = 0
      echo "Forcing file ncrcat mpi_status = $mpi_status"
   endif

   if ($mpi_status == 0 && $gr_stat != 0) then
      rm mycmdfile *.eo
   else
      echo 'ERROR in repackaging forcing (cpl history) files: See h*.eo, cmds_template, mycmdfile'
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

echo "------------------------"
if ($do_restarts == true) then
   echo "Restarts starts at "`date`
   
   cd rest

#    set pre_clean = false
#    set files_to_save = true
   set pre_clean = true
   set files_to_save = false

   foreach rd (`ls -d ${yr_mo}-*`)
      # Purge restart directories which don't fit in the frequency 
      # defined in setup_* by save_every_Mth_day_restart.
#       echo ' '
#       echo Contents of rest at the start of the dates loop
#       ls -dlt *

      echo ' '
      echo Processing $rd
      # Prevent archiving stray files.
      if (-f $rd) continue

      echo $rd | grep '[a-zA-Z]' 
      if ($status == 0) continue

      echo $rd | grep '\-'
      if ($status == 0) then
         set rd_date_parts = `echo $rd | sed -e "s#-# #g"`
         # set year        = $rd_date_parts[1]
         # set month       = $rd_date_parts[2]
         set day_o_month = $rd_date_parts[3]
         set sec_o_day   = $rd_date_parts[4]
         echo year = $year
      else
         echo continuing to next restart candidate.
         continue
      endif

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
         # set purge_date = $rd_date_parts[1]-$rd_date_parts[2]-$rd_date_parts[3]
         set purge_date = ${yr_mo}-${day_o_month}
         set weekday = `date --date="$purge_date" +%A`
         if ($weekday == $save_rest_freq) set purge = 'false'
   
      # Numbers can be tested inside the 'if' statement.
      else if (`echo $save_rest_freq | grep '[0-9]'`) then
         if ($day_o_month % $save_rest_freq == 0 && $sec_o_day == '00000') set purge = 'false'

      endif
   
      if ($purge == 'true') then
         echo "Removing restart directory $rd because it doesn't match "
         echo "         save_every_Mth_day_restart = $save_rest_freq and/or $sec_o_day != 00000"
# Moved into the purge script.
#>>>          rm -rf $rd
# ?      This generated an error message, but seems to have succeeded.
#        The message doesn't appear when this is done interactively on casper.
#        slurm.$jobid: "rmdir: failed to remove ‘2010-07-03-43200’: No such file or directory"
#        Also, those extra characters don't appear in interactive mode.
# ?      And why is the error message from rmdir, when I'm using 'rm -rf'?


      # This prevents the most recent (uncompressed) restart set from being archived
      # But this may be prevented more reliably by the selection in foreach ( $yr_mo).
      else 
         echo "Exporting restart file set ${rd} "
         echo "   to ${campaign}/${CASE}/rest/${year} "

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
    
         if ($status == 0) then
            grep tar *.eo | grep -v log >& /dev/null
            # grep failure = tar success = "not 0"
            set gr_stat = $status
         else
            # No eo files = failure of something besides tar.
            set gr_stat = 0
            echo "Restart file set, tar mpi_status = $mpi_status"
         endif
      
         if ($mpi_status == 0 && $gr_stat != 0) then
            rm tar*.eo
         else
            echo 'ERROR in repackaging restart files: See tar*.eo, mycmdfile'
            echo '      grep tar *.eo  yielded status '$gr_stat
            ls -l *.eo
            exit 20
         endif

      endif

      # Copy all the restart tars to campaign storage
      # 2019-4-26; Can't ever do from cheyenne batch nodes.
      if ($files_to_save == true) then
         # Remove the empty directory to prevent mv_to_campaign.csh from archiving it.
         rmdir $rd
   #       rm mycmdfile
   
         ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/rest/$yr_mo \
                                        ${campaign}/${CASE}/rest
         echo "WARNING: mv_to_campaign.csh ONLY SUBMITS THE REQUEST to globus"
         echo "         The directory must be manually removed after globus completes the copy."
      endif
   end

   #  case  time_str         local_dir                destination

   # OR         Print the command so that I can run it interactively after this job is done.
   # echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/rest/$year " \
   #      " ${campaign}/${CASE}/rest/$year"
   cd ..
   # Yes, only one level up for the rest directory.
endif


# 3) Obs space diagnostics.
#    Generated separately using ./obs_diags.csh

echo "------------------------"
if ($do_obs_space == true) then
   cd esp/hist
   echo " "
   echo "Location for obs space is `pwd`"

   # A csh pattern that will find the days of desired obs_seq files.
   # set obs_times_pattern = '2017-{01}-01'
   # Non-standard set;
   # set obs_times_pattern = '2017-{01}-{0*,1[0-4]}'
   set obs_times_pattern = $yr_mo

   # Create the obs lists and the Obs_seq directory, which diags_batch.csh needs to use.
   # Decompress, if needed.
   # For now, do this serially, since we've greatly reduced the size of the obs_seq files,
   # and won't compress them during the assimilation, so this won't be needed often.
   # Later I may want to expand compress.csh to handle this set of decompression (a date span).
# ? TURN THIS ON FOR PRODUCTION RUN
#    echo Making raw.list
#    ls -1 ${CASE}.dart.e.cam_obs_seq_final.${obs_times_pattern}* >! raw.list
# 
#    if (-f Obs_seqs.list) rm Obs_seqs.list
#    touch Obs_seqs.list
#    
#    echo Filling Obs_seqs.list
#    foreach f (`cat raw.list`)
#       if ($f:e == 'gz') then
#          gunzip $f
#          echo "../$f:r" >> Obs_seqs.list
#       else 
#          echo "../$f"   >> Obs_seqs.list
#       endif
#    end
# 
#    if (-d Obs_seqs) rm -rf Obs_seqs
#    echo Making Obs_seqs
#    mkdir Obs_seqs
#    cd Obs_seqs
#    ln -s ../${CASE}.dart.e.cam_obs_seq_final.${obs_times_pattern}* .
#    ls -l
#    cd ..
#    
   # Harvest dates for the diagnostics directories used by diags_batch. 
   # set ofile = `head -n 1 Obs_seqs.list `
   # set first = `echo $ofile:e | sed -e "s#-# #g"`
   # set ofile = `tail -n 1 Obs_seqs.list `
   # set last  = `echo $ofile:e | sed -e "s#-# #g"`
   # set date_span = $first[1].$first[2].$first[3]-$last[1].$last[2].$last[3]
   
#    ${CASEROOT}/diags_batch.csh ${CASE} Obs_seqs $obs_space
#    if ($status == 0) then
#       # Done with obs_seq_final files.  Prep them for archiving.
         tar -z -c -f ${CASE}.cam_obs_seq_final.${yr_mo}.tgz \
               ${CASE}.dart.e.cam_obs_seq_final.${obs_times_pattern}* &
#    else
#       echo diags_batch.csh failed with status = $status
#       exit 35
#    endif
# 
# ? TURN THIS ON FOR PRODUCTION RUN
# # Generate pictures of the obs space diagnostics using matlab.
#    cd $obs_space
# 
#    echo "addpath('$DART/diagnostics/matlab','-BEGIN')" >! script.m
#    echo "fname = 'obs_diag_output.nc';"                >> script.m
# 
#    foreach copy (totalspread bias)
#    foreach func (plot_rmse_xxx_evolution plot_rmse_xxx_profile)
#       echo "copystring = '$copy';"           >> script.m
#       echo "$func(fname,copystring)"         >> script.m
#    end
#    end
# 
#    matlab -r script
#    
#    if ($status == 0) then
#       # Prep the obs space diagnostics for archiving;
#       # obs_diag_output.nc and matlab output.
#       cd ..
#       tar -z -c -f ${obs_space}.tgz $obs_space 
#    else
#       echo ERROR: matlab failed with status = $status.
#       # Let the tar of obs_seq files finish.
#       wait
#       exit
#    endif
#    
   # Move the obs space diagnostics to $project.
   set obs_proj_dir = ${project}/${CASE}/esp/hist
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
   if (  -f ${CASE}.cam_obs_seq_final.${obs_times_pattern}.tgz && \
       ! -z ${CASE}.cam_obs_seq_final.${obs_times_pattern}.tgz) then
      mv ${CASE}.cam_obs_seq_final.${obs_times_pattern}.tgz $obs_proj_dir
      rm ${CASE}.dart.e.cam_obs_seq_final.${obs_times_pattern}*
      rm -rf Obs_seqs
      echo "Moved ${CASE}.cam_obs_seq_final.${obs_times_pattern}.tgz "
      echo "   to $obs_proj_dir"
   else
      echo "${CASE}.cam_obs_seq_final.${obs_times_pattern}.tgz cannot be moved"
      exit 35
   endif

   cd ../..
   
endif

#--------------------------------------------

# 4) CESM history files
#    Needs to come before 5) State space, so that tarred h# files
#      are created in $project and h# files are deleted
#      before the state space mv_to_campaign copies everything in atm/hist 
#      to campaign storage.
#    + CLM output (.h1.) (for Lombardozzi) may need monthly averaging.
#    + .h0. does not., but she only needs daily; concatenate the -00000 files.
#    Should all members be saved?  Then use cmdfile.
#    Or just 1?  Much simpler
#    See earlier sections for generating a list of files.

echo "------------------------"
if ($do_history == true) then
   # If cam files are too big, do them separately and monthly.
   set m = 1
   echo "There are $#components components (models)"
   while ($m <= $#components)
      ls $components[$m]/hist/*h0* >& /dev/null
      if ($status == 0) then
         cd $components[$m]/hist
         echo " "
         echo "Location for history is `pwd`"
      else
         echo "Skipping $components[$m]/hist"
         @ m++
         continue
      endif

      set i = 1
      while ($i <= $ensemble_size)
         set NINST = `printf %04d $i`
         set inst_dir = ${project}/${CASE}/$components[$m]/hist/${NINST}
         set yearly_file = ${CASE}.$models[$m]_${NINST}.TYPE.${year}.nc

         if (! -d $inst_dir) then
            mkdir -p $inst_dir
         else if (-f ${inst_dir}/$yearly_file) then
            mv ${inst_dir}/$yearly_file . 
         else
            echo "ERROR: $inst_dir exists, but does not have $yearly_file in it"
            exit 50
         endif

         @ i++
      end
   
      # Make a cmd file to append this month's time series to the yearly file in $project
      # Start with a template of all the instances of one file type.
   
      # This is what's in cmds_template, which is re-used here.
      #       set inst_dir = ${project}/${CASE}/cpl/hist/${NINST}
      #       set yearly_file = ${CASE}.cpl_${NINST}.TYPE.${year}.nc
      #       echo "mv ${inst_dir}/$yearly_file . ; ncrcat --rec_apn    $yearly_file " \
      #            "${CASE}.cpl_${NINST}.TYPE.${yr_mo}-*.nc ${inst_dir}/$yearly_file &> " \
      #            "TYPE_${NINST}.eo " \
      set cmds_template = cmds_template_$models[$m]
      if (-f ../../cpl/hist/cmds_template) then
         sed -e "s#cpl_#$models[$m]_#g;s#cpl#$components[$m]#"  \
             ../../cpl/hist/cmds_template >! $cmds_template
      else
         echo "ERROR: ../../cpl/hist/cmds_template is missing; need it for archiving h# files."
         echo "       It should have been created in section 1 of this script."
         exit 50
      endif

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
         set dates = `ls *0001.h${type}.${yr_mo}-*`
         if ($#dates == 0) break
   
         # There are ensemble_size commands in cmds_template.
         if ($models[$m] == 'cam' && $type == 0) then
            # If cam.h0 ends up with more than PHIS, comment this out
            # and fix the h0 purging in the state_space section.
            # Actually; skip cam*.h0. because of the purging done by assimilate.csh.
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
         echo "Skipping $models[$m]/hist because $mycmdfile has size 0"
      else
         echo "   history mpirun launch_cf.sh starts at "`date`
         mpirun -n $tasks ${CASEROOT}/launch_cf.sh ./${mycmdfile}
         set mpi_status = $status
         echo "   history mpirun launch_cf.sh ends at "`date`
      
         ls *.eo > /dev/null
         if ($status == 0) then
            grep ncrcat *.eo >& /dev/null
            # grep failure = ncrcat success = "not 0"
            set gr_stat = $status
         else
            echo "No eo files = failure of something besides g(un)zip."
            echo "   History file ncrcat mpi_status = $mpi_status"
            set gr_stat = 0
         endif
      
         if ($mpi_status == 0 && $gr_stat != 0) then
            rm ${cmds_template}  ${mycmdfile} *.eo
         else
            echo "ERROR in repackaging history files: See $components[$m]/hist/"\
                 'h*.eo, cmds_template*, mycmdfile*'
            echo '      grep ncrcat *.eo  yielded status '$gr_stat
            exit 20
         endif
      endif

      cd ../..

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
   cd esp/hist
   echo " "
   echo "Location for state space is `pwd`"
   
   mkdir $yr_mo
   foreach stage ($stages_all)
      set ext = e
      if ($stage == output) set ext = i
   
      # Ignoring posterior inflation for now
      foreach stat  (mean sd priorinf_mean priorinf_sd)
         echo $stat | grep inf 
         if ($status == 0) set ext = rh
         echo "stage, stat, ext = $stage $stat $ext"
   
         # Do not want all stages of inflation. 
# >>> ?    Except that the user can set save_all_inf = TRUE in assimilate.csh
#          if (($stat == 'priorinf_mean' || $stat == 'priorinf_sd') && \
#              $stage == 'preassim') then
#             rm ${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}-*
#          else 
            if (! -f     ${yr_mo}/${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}.tar) then
               tar -c -f ${yr_mo}/${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}.tar  \
                                  ${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}-*
               if ($status == 0) then
                  rm ${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}-*
               else
                  echo tar -c -f ${CASE}.dart.${ext}.cam_${stage}_${stat}.${yr_mo}.tar failed
                  exit 30
               endif
            endif
#          endif
      end
   end
   # Echoing the command instead of running it, so that this script can run on cheyenne
   # instead of casper.
   # FIX as below before using this
   #    echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/esp " \
   #         " ${campaign}/${CASE}/esp/hist"

   ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/esp/hist/${yr_mo}  \
                                  ${campaign}/${CASE}/esp/hist
 
 #    + atm/hist/$inst.e.$stages_except_output
 #      compressed already  73 Mb * 80 * 5 dates =  25 Gb/mo.
 #      compress?  85 Mb * 80 * 120 dates = 816 Gb -> 710 Gb/mo.
 #                              ^ f*00[1-4] saved all (120/mo) times, due to missing .e.
 #                                from save_stages_freq archiving section.
 #      Already done during assim.  The savings is over a Tb/year.
   cd ../../atm/hist
   
   # Archive the members, which are here and not in esp/hist (where the means are)
   set files = `ls ${CASE}.cam_0001.e.$stages_all[1].${yr_mo}*`
   echo "Files from which atm preassim allinst dates will be gathered:"
   echo "  $files"
   if ($#files == 0) then
      echo "There are no .e.$stages_all[1].${yr_mo} files in atm/hist.  Continuing."
   else

      mkdir $yr_mo
      set dates = ()
      foreach f ($files)
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
#>>>
#             echo "For now, do not rm ${CASE}.cam_[0-9]*.e.${stage}.${date}*"
            # Preserve one of the CAM h0 files.
            # This is apparently very tricky.  The interaction of the assimilate.csh purging 
            # of h0 files with globus leads to globus errors (file not found).
            # mv    ${CASE}.cam_0001.h0.${date}* ..
            # rm    ${CASE}.cam_*.h0.${date}*
            # mv ../${CASE}.cam_0001.h0.${date}*  .
         else
            echo tar -c -f ${CASE}.cam_allinst.e.$stage.${date}.tar failed
            exit 40
         endif
      endif
      end
      end
      
    # Echoing the command instead of running it, so that this script can run on cheyenne
    # instead of casper.
    # FIX as below before using this
    # echo "${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/atm " \
    #      " ${campaign}/${CASE}/atm/hist"

      ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/atm/hist/$yr_mo ${campaign}/${CASE}/atm/hist
   
   endif

   # Archive DART log files (and others?)

   cd ../../logs
   
   # Create a list of files to archive.
   # Don't use cmdfile; to prevent confusion in the list.
   set list = ()
   foreach f (`ls da.log*`)
      if ($f:e == 'gz') then
         gunzip $f 
         set f = $f:r
      endif
      grep -l "valid time of model is $year $month" $f > /dev/null
      set gstatus = $status
      if ($gstatus == 0) then
#          echo "   $f added to list"
         set list = ($list $f)
#       else
#          echo "   $f not added to list with status $gstatus"
      endif
   end
   echo "Archiving "
   echo $list

   # Logs from more components could be added here.
   if ($#list > 0) then
      if (! -d $yr_mo) mkdir $yr_mo
      tar -z -c -f ${yr_mo}/da.log.${yr_mo}.tar $list
      if ($status == 0) then
         rm $list
         ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${local_arch}/logs/${yr_mo} ${campaign}/${CASE}/logs
      else
         echo "Tar of da.logs of $yr_mo failed.  Not archiving them"
      endif
   else
      echo "da.log list has none in it"
   endif

   cd ../..

endif 

#--------------------------------------------

exit
#==========================================================================
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
                                                              
