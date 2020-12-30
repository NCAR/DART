#!/bin/tcsh -f

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Script to set up the submission of a variety of jobs:
#    full cycle(s)
#    assim only
#    flexible numbers of cycles/job and resubmissions
# and to check whether the case is ready for that submission:
#    prevent running assimilations when there are too many cesm.log files in $rundir
#    check the requested start date against the one in rpointer.atm_0001
#
# After all the desired changes have been applied, preview_namelists will be run
# as this is necessary to generate the namelists that will be used for the next 
# CESM advance. A little post-processing is done on Buildconf/[cam,clm].input_data_list
# to generate new files that are sorted and have unique entries. That file is what should
# be under version control.
#
# Run in $CASEROOT

if ($#argv != 4) then
   echo "Usage: submit first_date last_date cycles_per_job queue"
   echo "       first_date, last_date = YYYY-MM-DD-SSSSS"
   echo "       first_date = last_date -> assimilation only."
   echo "       cycles_per_job: The number of jobs will be calculated from the dates"
   echo "                       Set to 1 for assimilation only jobs."
   echo "       If 'queue' is omitted, the wall clock will be calculated, then exit"
   exit 10
endif

# Get CASE environment variables from the central variables file.
source ./data_scripts.csh
echo "data_CASEROOT    = ${data_CASEROOT}"
echo "data_NINST       = ${data_NINST}"
echo "data_scratch     = ${data_scratch}"
echo "data_CESM_python = ${data_CESM_python}"

set first_date     = $1
set last_date      = $2
set cycles_per_job = $3
set queue          = $4

# Check whether there are too many cesm.log files in rundir.
set num_logs = `ls -1 ${data_scratch}/run/cesm.log* | wc -l`

if ($first_date == $last_date && $num_logs > 3 || \
    $first_date != $last_date && $num_logs > 2 ) then
   echo "ERROR: too many cesm.log files in ${data_scratch}/run; remove the extraneous"
   exit 20
endif

# Make sure that requested initial date is what CAM will use.
if (-f ${data_scratch}/run/rpointer.atm_0001) then
   grep $first_date ${data_scratch}/run/rpointer.atm_0001
   if ($status != 0) then
#     $first_date != $last_date && 
# ?   This relies on stage_cesm_files.template, which is currently NOT made by setup_advanced.
      sed -e "s#NO-DATE-YET#$first_date#" stage_cesm_files.template >! stage_cesm_files
      if ($status != 0) then
         echo "ERROR: rpointer.atm_0001 has the wrong date, but creating stage_cesm_files failed"
         exit 30
      endif
      echo "Running stage_cesm_files for $first_date"
      chmod 755 stage_cesm_files
      ./stage_cesm_files
      if ($status != 0) then
         echo "ERROR: stage_cesm_files failed."
         echo "       Check: restart_date, DOUT_S, CONTINUE_RUN, ..."
         exit 40
   endif
endif

# Translate dates into cycles.
set start = `echo $first_date | sed -e "s#-# #g"`
set end   = `echo $last_date  | sed -e "s#-# #g"`

set days_in_mo = (31 28 31 30 31 30 31 31 30 31 30 31)
if (($start[1] % 4) == 0 ) set days_in_mo[2] = 29

@ years  = $end[1] - $start[1]
@ months = $end[2] - $start[2]
@ days   = $end[3] - $start[3]
@ secs   = $end[4] - $start[4]
if ($secs < 0 ) then
   @ secs = $secs + 86400
   @ days--
endif
if ($days < 0 ) then
   @ days = $days + $days_in_mo[$start[2]]
   @ months--
endif
if ($months < 0 ) then
   @ months = $months + 12
   @ years--
endif

if (($end[4] < $start[4] && $end[4] != '00000') || \
    ($end[4] > $start[4] && $end[3] != $start[3]) ) then
   echo "ERROR: submit requires an integral number of days "
   echo "       or starting and ending on the same day"
   echo "       or ending on 00000 of a future day"
   exit 50
endif

# Make sure docn will exit if the hindcast span is not in the data file.
set sst_use = `grep taxmode user_nl_docn_0001`
set sst_mode = `echo $sst_use[3] | sed -e 's#"##g'`
# echo "$sst_use"
# echo "sst_mode is $sst_mode"
if ($sst_mode != limit) then
   echo 'In user_nl_docn* taxmode should be "limit"'
   exit 60
endif


# Turn off short term archiver for the new month,
# after it was activated for the last job of a month.
# Change; now that we're trying to run whole months,
# it's convenient to make st_archive run automatically
# at the end of successful jobs.
# if ($start[3] == "01") ./xmlchange DOUT_S=FALSE
if ($end[3] == "01" && $end[4] == "00000") then
   ./xmlchange DOUT_S=TRUE
else
   ./xmlchange DOUT_S=FALSE
endif

# I won't run more than a month, so I'm leaving the years
# contribution (0) out of this tally.
@ cycles = ( $days * 4 ) + ( $secs / 21600 )
if ($months > 0) @ cycles = $cycles + ( $days_in_mo[$start[2]] * 4 )

@ resubmissions = ( $cycles / $cycles_per_job ) - 1
if ($resubmissions < 0 && $first_date != $last_date ) then
   echo "ERROR: resubmissions < 0, so cycles_per_job is probably too big"
   exit 70
endif
./xmlchange DATA_ASSIMILATION_CYCLES=${cycles_per_job},RESUBMIT=${resubmissions}
echo "$cycles cycles will be distributed among $resubmissions +1 jobs"

# 10 minutes is what's needed for the first cycle,
# even if cycles = 0 for an assimilation only job,
# but later cycles need more.  It doubles in ~60 cycles.
# >>> This may be fixed with Brian Dobbins > remove the nonlinear term.
# @ job_minutes = 10 + ( $cycles_per_job * ( 10 + (( $cycles_per_job * 10) / 70 )))
# 12 hours: @ job_minutes = 720

@ job_minutes = ( $cycles_per_job * 7 ) + 10 
# @ job_minutes = 12 * 60
if ($job_minutes > 720) set job_minutes = 720
@ wall_hours  = $job_minutes / 60
@ wall_mins   = $job_minutes % 60
set wall_time = `printf %02d:%02d $wall_hours $wall_mins`
echo "Changing run time to $wall_time in env_batch.xml"

if ($#argv == 3) then
   echo "Seeing if time span will fit in 12:00:00"
   exit 80
endif
if ($job_minutes > 960) then
   echo "ERROR: too many cycles requested.  Limit wall clock is 12:00:00"
   exit 90
endif

./xmlchange --subgroup case.run --id JOB_WALLCLOCK_TIME      --val ${wall_time}:00
./xmlchange --subgroup case.run --id USER_REQUESTED_WALLTIME --val ${wall_time}
./xmlchange --subgroup case.run --id JOB_QUEUE            --val $queue
./xmlchange --subgroup case.run --id USER_REQUESTED_QUEUE --val $queue

# Choose a version of case_run.py to use; with(out) a CAM run.
cd ${data_CESM_python}/case
if (-l case_run.py) then
   rm case_run.py
else
   echo 'ERROR: case_run.py is not a link.  Make it one'
   exit 100
endif

echo "Comparing dates for linking"
if ("$first_date" == "$last_date" ) then
   # Assimilation only; link an assim-only copy to the expected name
   ln -s case_run_only_assim.py case_run.py
   set init_files = `wc -l ${data_scratch}/run/cam_init_files`
   echo "Checking numbers of files " $init_files[1] ${data_NINST}
   if ( $init_files[1] != ${data_NINST} ) then
      echo "ERROR: the forecast didn't finish; not enough files in cam_init_files"
      exit 110
   endif
else
   # Hindcast + assimilation; link the right version to the expected name.
   ln -s case_run_cam+assim.py case_run.py
endif   
echo "case_run.py is linked to"
ls -l case_run.py
cd -

# Generate the namelists that will be used for the model advance and should
# be version-controlled. Also generates metadata that needs post-processing
# before being version-controlled.
# The *modelio* namelists will get regenerated as a matter of course
# because preview_namelists is always run by the first cycle.
# The *modelio* namelists have logfile names with unique strings that
# are generated when preview_namelists runs. Consequently, the logfile
# names generated here will be changed when preview_namelists is run by
# the first cycle. The --skip-preview-namelist option to case.submit prevents
# running preview_namelists for each of the DATA_ASSIMILATION_CYCLES,
# which is a colossal waste of time.

./preview_namelists

echo "Creating Buildconf/[cam,clm].input_data_list.sorted files."

cd Buildconf
sort cam.input_data_list | uniq > cam.input_data_list.sorted
sort clm.input_data_list | uniq > clm.input_data_list.sorted
cd -

echo "If this is the start of a new month, issue a pull request to DART_CASES:"
echo '   % git status -uno | sed -e "s#modified:#git add#" >! push_prep.csh'
echo '   Edit push_prep.csh to make it a csh script '
echo '   % csh push_prep.csh'
echo '   % git commit '
echo "     with comments about the important modifications."
echo '   % git push origin '$data_CASEROOT:t
echo "   On github...kdraeder/DART_CASES issue the pull request."
# Echo the submit command, without generating new ${comp}_in_#### files.
echo "After approval, submit the job using"
echo './case.submit -M begin,end --skip-preview-namelist'

exit 0

