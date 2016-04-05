#/bin/csh
#
# This code may (or may not) be part of the NOGAPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
# Function:
#     script to stage a DART/NOGAPS assimilation experiment 
#
#############################################################################

# Set the base directories
# FIXME: scratch_dir sounds too temporary for something that's really
# a working dir and needs to persist between invocations of filter.

# set dart_base_dir     = "/home/coral/hansenj/DART_experiments/cookbook/DART/models/NOGAPS"
set dart_base_dir     = "/fs/image/home/thoar/SVN/DART/models/NOGAPS"
set ocards_files      = "${dart_base_dir}/templates/ocards_files"
set climo_dir         = "${dart_base_dir}"
set scratch_dir       = "/ptmp/thoar/NOGAPS1"
set experiment_name   = "test1"
set experiment_dir    = ${scratch_dir}/${experiment_name}

# Set the default date.

set dtg = 2008080100

# NOGAPS executable

set resolution        = 159 
set n_levels          = 30 
#set NOGAPS_exec_dir   = "/home/coral/hansenj/nogaps_T${resolution}L${n_levels}"
set NOGAPS_exec_dir   = "${dart_base_dir}/src"
set NOGAPS_exec_name  = "got${resolution}l${n_levels}"
set perturb_dir       = "/home/coral/hansenj/DART_experiments/cookbook/DART/models/NOGAPS/init_perts_T${resolution}L${n_levels}"
set climo             = "${scratch_dir}/climo${resolution}"

# archive directories  

set archive              = "newton"
set archive_ens_dir      = "./outp"
set archive_climo_dir    = "./"
set archive_analysis_dir = "./analyses_T${resolution}L${n_levels}"

# DART executables and scripts

set DART_exec_dir          = "${dart_base_dir}/work"
set DART_script_dir        = "${dart_base_dir}/shell_scripts"

set DART2NOGAPS_exec_name  = "${DART_exec_dir}/dart_to_nogaps"
set NOGAPS2DART_exec_name  = "${DART_exec_dir}/nogaps_to_dart"
set PERFECT_exec_name      = "${DART_exec_dir}/perfect_model_obs"
set FILTER_exec_name       = "${DART_exec_dir}/filter"
set WAKEUPFILTER_exec_name = "${DART_exec_dir}/wakeup_filter"
set TRANSTIME_exec_name    = "${DART_exec_dir}/trans_time"
set OBSTOOL_exec_name      = "${DART_exec_dir}/obs_sequence_tool"
set ADVTIME_exec_name      = "${DART_exec_dir}/advance_time"

# observation files directory

set DART_obs              = "/home/coral/hansenj/obs_data"


# NOGAPS log-file directory

set NOGAPS_log_dir = "${scratch_dir}/log"


# ending lead time and wall clock 

set endtau = '1'          # this is overwritten by the actual advance time
set clock_time = '0:30'


# various MPI parameters - It is important that the total jproc*iproc
# matches the number of MPI tasks you run with, even for the converters.

set iproc             = 1

# NOTE: The original idea was to have the runscripts be general
#       and set anything specific in this config script.  However,
#       at the moment the runscripts are specific so the MPI 
#       variable isn't used.
#
# MPI   must be set to the exact command to launch an MPI job -
#       a trailing space is good protection.

if ( $?LSB_QUEUE ) then
   # We must be using an LSF-based system.
   # count the number of hosts, which equals the number of MPI tasks.
   # jproc must match that.  for a batch system other than LSF either
   # add a section here, or set jproc by hand.
   set jproc = `echo $LSB_HOSTS | wc -w`
   set   MPI = 'mpirun.lsf '
else if ( $?PBS_QUEUE ) then
   # We must be using a PBS-based system.
   # count the number of hosts, which equals the number of MPI tasks.
   # jproc must match that.  for a batch system other than LSF either
   # add a section here, or set jproc by hand.
   set jproc = `echo $PBS_HOSTS | wc -w` # TJH FIXME
   set   MPI = 'mpirun '
else
   set jproc = 8 
   set   MPI = 'chokeme '
endif

set jsplit=2
@ nproc = $iproc * $jproc

#############################################################################
# print a summary table if you like
#############################################################################

if ( 1 == 2 ) then
   echo "ADVTIME_exec_name is $ADVTIME_exec_name"
   echo "archive is $archive"
   echo "archive_analysis_dir is $archive_analysis_dir"
   echo "archive_climo_dir is $archive_climo_dir"
   echo "archive_ens_dir is $archive_ens_dir"
   echo "climo_dir is $climo_dir"
   echo "clock_time is $clock_time"
   echo "DART2NOGAPS_exec_name is $DART2NOGAPS_exec_name"
   echo "dart_base_dir is $dart_base_dir"
   echo "DART_exec_dir is $DART_exec_dir"
   echo "DART_obs is $DART_obs"
   echo "DART_script_dir is $DART_script_dir"
   echo "dtg is $dtg"
   echo "endtau is $endtau"
   echo "experiment_dir is $experiment_dir"
   echo "experiment_name is $experiment_name"
   echo "FILTER_exec_name is $FILTER_exec_name"
   echo "iproc is $iproc"
   echo "jproc is $jproc"
   echo "jsplit is $jsplit"
   echo "MPI is $MPI"
   echo "n_levels is $n_levels"
   echo "NOGAPS2DART_exec_name is $NOGAPS2DART_exec_name"
   echo "NOGAPS_exec_dir is $NOGAPS_exec_dir"
   echo "NOGAPS_exec_name is $NOGAPS_exec_name"
   echo "NOGAPS_log_dir is $NOGAPS_log_dir"
   echo "OBSTOOL_exec_name is $OBSTOOL_exec_name"
   echo "ocards_files is $ocards_files"
   echo "PERFECT_exec_name is $PERFECT_exec_name"
   echo "perturb_dir is $perturb_dir"
   echo "resolution is $resolution"
   echo "scratch_dir is $scratch_dir"
   echo "TRANSTIME_exec_name is $TRANSTIME_exec_name"
   echo "WAKEUPFILTER_exec_name is $WAKEUPFILTER_exec_name"
endif
   
exit 0
   
