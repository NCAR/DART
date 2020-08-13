#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Script to assimilate observations using DART and the POP ocean model.
# This presumes two directories exists that contain all the required bits
# for POP and for DART.
#
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluefire'
#
# the normal way to submit to the queue is:    bsub < run_filter
#
# an explanation of the most common directives follows:
# -J Job name (master script job.csh presumes filter_server.xxxx.log)
# -o STDOUT filename
# -e STDERR filename
# -P      account
# -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
# -n number of processors  (really)
##=============================================================================
#
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -q regular
#BSUB -n 16
#BSUB -R "span[ptile=2]"
#BSUB -P 86850054
#BSUB -W 2:00
#BSUB -N -u ${USER}@ucar.edu
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
## PBS is used on the CGD Linux cluster 'bangkok'
## PBS is used on the CGD Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub run_filter
##
## an explanation of the most common directives follows:
## -N     Job name
## -r n   Declare job non-rerunable
## -e <arg>  filename for standard error
## -o <arg>  filename for standard out
## -q <arg>   Queue name (small, medium, long, verylong)
## -l nodes=xx:ppn=2   requests BOTH processors on the node. On both bangkok
##                     and calgary, there is no way to 'share' the processors
##                     on the node with another job, so you might as well use
##                     them both. (ppn == Processors Per Node)
##=============================================================================
#
#PBS -N filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.log
#PBS -q medium
#PBS -l nodes=8:ppn=2

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv MPI         mpirun.lsf

else if ($?PBS_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by PBS
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST
   setenv MPI         mpirun

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     POP
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv MPI         csh

endif

#----------------------------------------------------------------------
# Just an echo of the job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $HOST"
echo "${JOBNAME} ($JOBID) started      at "`date`
echo

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#----------------------------------------------------------------------

setenv TMPDIR /ptmp/${user}/${JOBNAME}/job_${JOBID}

mkdir -p ${TMPDIR}
cd ${TMPDIR}

set CENTRALDIR = `pwd`
set myname = $0          # this is the name of this script

# some systems don't like the -v option to any of the following 

set OSTYPE = `uname -s`
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -vp'
      setenv   MOVE 'mv -fv'
      breaksw
endsw

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = /fs/image/home/${user}/SVN/DART/models/POP
set  POPDIR = /ptmp/${user}/POP/osse
set  OBSERVATIONDIR = /ptmp/${user}/POP_OSSE/1_7_Jan2000

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
#-----------------------------------------------------------------------------

# executables
 ${COPY} ${DARTDIR}/work/filter                     .
 ${COPY} ${DARTDIR}/work/wakeup_filter              .
 ${COPY} ${DARTDIR}/work/dart_to_pop                .
 ${COPY} ${DARTDIR}/work/pop_to_dart                .
 ${COPY} ${DARTDIR}/work/restart_file_tool          .

# shell scripts
 ${COPY} ${DARTDIR}/shell_scripts/cesm1_x/advance_model.csh .

# data files
 ${COPY} ${DARTDIR}/work/input.nml                  .
 ${COPY} ${OBSERVATIONDIR}/obs_seq.out              .

#-----------------------------------------------------------------------------
# Get the POP executable, control files, and data files.
# trying to use the CCSM naming conventions
#-----------------------------------------------------------------------------

 ${COPY} ${POPDIR}/pop                       .
 ${COPY} ${POPDIR}/pop_in.part1              .
 ${COPY} ${POPDIR}/pop_in.part2              .

 ${COPY} ${POPDIR}/gx3v5_tavg_contents       .
 ${COPY} ${POPDIR}/gx3v5_movie_contents      .
 ${COPY} ${POPDIR}/gx3v5_history_contents    .
 ${COPY} ${POPDIR}/gx3v5_transport_contents  .

 ${COPY} ${POPDIR}/vert_grid.gx3v5              .
 ${COPY} ${POPDIR}/horiz_grid.gx3v5.r8ieee.le   .
 ${COPY} ${POPDIR}/topography.gx3v5.i4ieee.le   .

#-----------------------------------------------------------------------------
# Determine the number of ensemble members from input.nml,
# It may exist in more than one place - we will use the first instance.
# Parse out the filter_nml string and use the next hunk of lines.
# ditto for the advance command
#-----------------------------------------------------------------------------

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set  ADVANCESTRING = `grep -A 42 filter_nml input.nml | grep adv_ens_command`
set  ensemble_size = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`
set        ADV_CMD = `echo  $ADVANCESTRING[3] | sed -e 's#,##' -e 's#"##g'`

echo "The model advance command is ${ADV_CMD}"

#-----------------------------------------------------------------------------
# detect whether the model is supposed to run as an MPI job or not
# by reading the "async = " from the &filter_nml namelist in input.nml.
#-----------------------------------------------------------------------------

set ASYNCSTRING = `grep -A 42 filter_nml input.nml | grep async`
set  ASYNC_TYPE = `echo $ASYNCSTRING[3] | sed -e 's#,##'`

if ( "${ASYNC_TYPE}" == "0" || "${ASYNC_TYPE}" == "2") then
  set parallel_model = "false"
  echo "The model is believed to be single-threaded."
else if ( "${ASYNC_TYPE}" == "4") then
  set parallel_model = "true"
  echo "The model is believed to be MPI-aware."
else
  echo 'ERROR - Cannot autodetect async value in the filter_nml namelist in input.nml.'
  echo 'ERROR - hardcode the parallel_model shell variable and comment out these lines.'
  exit -1
  set parallel_model = "false"
endif

#-----------------------------------------------------------------------------
# Block 1: convert N POP restart files to DART initial conditions file(s).
# Since the initial ensemble may not all have the desired timestamp, we
# will use restart_file_tool to use a consistent date in the header of
# all the DART initial conditions files. At the end of this block, 
# we have DART restart files   filter_ics.[1-N]    that
# came from pointer files    rpointer.ocn.[1-N].restart
#
# DART requires that POP uses pointer files and that the POP restart files
# are netCDF format. The experiment should be initialized such that there
# are "ensemble_size" number of POP restart files and matching pointer files.
# The pointer files should have the absolute path to the restart file.
#
# DART namelist settings appropriate/required:
# &filter_nml:           restart_in_file_name    = 'filter_ics'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &pop_to_dart_nml:      pop_to_dart_output_file = 'dart_ics',
#
# &restart_file_tool_nml: <see list that follows>
#  input_file_name              = "dart_input",
#  output_file_name             = "dart_output",
#  ens_size                     = 1,
#  single_restart_file_in       = .true.,
#  single_restart_file_out      = .true.,
#  overwrite_data_time          = .true.,
#  overwrite_advance_time       = .true.,
#  new_data_days                = 144731,       [1 january 2000]
#  new_data_secs                =      0,       [midnight]
#  input_is_model_advance_file  = .false.,
#  output_is_model_advance_file = .false.,
#  gregorian_cal                = .true.
#  new_advance_days             =  -1,
#  new_advance_secs             =  -1
#-----------------------------------------------------------------------------
# ensure namelists have desired values ...
#-----------------------------------------------------------------------------

# Gregorian 1 Jan 2000 <==> DART 145731

echo ':0'                               >! ex_commands
echo '/restart_file_tool_nml'           >> ex_commands
echo '/write_binary_restart_files'      >> ex_commands
echo ':s/.false./.true./'               >> ex_commands
echo '/overwrite_data_time'             >> ex_commands
echo ':s/.false./.true./'               >> ex_commands
echo '/new_data_days'                   >> ex_commands
echo ':s/-1/145731/'                    >> ex_commands
echo '/new_data_secs'                   >> ex_commands
echo ':s/-1/0/'                         >> ex_commands
echo ':wq'                              >> ex_commands

( ex input.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

cat pop_in.part1 pop_in.part2 >! pop_in

set member = 1
while ($member <= $ensemble_size)

   # grab the POP pointer file and dereference it 
   # Copy the POP restart file ... we will be updating it.
   # then link the POP restart file to the name for 'pop_to_dart'
   ${COPY} ${POPDIR}/rpointer.ocn.$member.restart .
   set OCN_RESTART_FILENAME = `head -1 rpointer.ocn.$member.restart`
   ${COPY} ${POPDIR}/${OCN_RESTART_FILENAME} .
   ln -sf ${OCN_RESTART_FILENAME} pop.r.nc

#  echo "Changing iyear of ${OCN_RESTART_FILENAME} to 2000"
#  ncatted -O -h -a iyear,global,o,l,2000 ${OCN_RESTART_FILENAME}

   ./pop_to_dart || exit 1

   ${MOVE} dart_ics dart_input

   ./restart_file_tool || exit 1

   # set the filename expected by DART for the initial conditions 
   set DART_IC_FILE = `printf filter_ics.%04d $member`
 
   ${MOVE} dart_output ${DART_IC_FILE}

   @ member++ 
end

#-----------------------------------------------------------------------------
# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LSB_QUEUE || $?PBS_QUEUE) then

    # Must be using LSF or PBS as the queueing system.
    echo "Using ${MPI} for execution"

    # each filter task advances the ensembles, each running on 1 proc.
    if ( "$parallel_model" == "false" ) then

       ${MPI} ./filter

    else

    # 1) filter runs in parallel until time to do a model advance.
    # 2) advance_model.csh successively runs N POP instances,
    #    each using the entire processor set.
    # 3) wakeup_filter wakes up filter so it can continue.

      \rm -f model_to_filter.lock filter_to_model.lock
      mkfifo model_to_filter.lock filter_to_model.lock

      set filterhome = ~/.filter$$
      if ( ! -e $filterhome) mkdir $filterhome

      # this starts filter but also returns control back to
      # this script immediately.

      ( setenv HOME $filterhome; ${MPI} ./filter ) &

      while ( -e filter_to_model.lock )

        set todo=`cat < filter_to_model.lock`
        echo "todo received, value = ${todo}"

        if ( "${todo}" == "finished" ) then
          echo "main script: filter done."
          wait
          break

        else if ( "${todo}" == "advance" ) then

          echo "calling model advance now:"
          ./advance_model.csh 0 ${ensemble_size} filter_control00000 || exit 9

          echo "restarting filter."
          ${MPI} ./wakeup_filter

        else

          echo "main script: unexpected value received."
          break

        endif

      end

      echo "filter finished, removing pipes."
      \rm -f model_to_filter.lock filter_to_model.lock

      if ( -d $filterhome) rmdir $filterhome
    endif

else

    # If you have a linux cluster with no queuing software, use this
    # section. The list of computational nodes is given to the mpirun
    # command and it assigns them as they appear in the file. In some
    # cases it seems to be necessary to wrap the command in a small
    # script that changes to the current directory before running.

    echo "running with no queueing system"
    echo "This is untested for POP -- ending now."
    exit

    # before running this script, do this once. the syntax is
    # node name : how many tasks you can run on it
    #setenv MYNODEFILE ~/nodelist
    #echo "node7:2" >! $MYNODEFILE
    #echo "node5:2" >> $MYNODEFILE
    #echo "node3:2" >> $MYNODEFILE
    #echo "node1:2" >> $MYNODEFILE

#   for compas
    setenv NUM_PROCS `cat nodelist-pgi | wc -l`
    set MPIRUN = /opt/mpich/myrinet/pgi/bin/mpirun
    setenv MPI "$MPIRUN -np $NUM_PROCS -nolocal -machinefile nodelist-pgi"

#   for atlas-pgi
    setenv NUM_PROCS `cat nodelist-pgi | wc -l`
    set MPIRUN = /share/apps/mpich1/pgi/bin/mpirun
    setenv MPI "$MPIRUN -np $NUM_PROCS -nolocal -machinefile nodelist-pgi"

#   for atlas-gfortran
    set MPIRUN = /share/apps/openmpi/gfortran/bin/mpirun
    setenv MPI "$MPIRUN --hostfile nodelist-gfortran --mca mtl mx --mca pml cm -np 72"

    echo "MPI = ${MPI}"

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the mitgcmuv jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

    \rm -f model_to_filter.lock filter_to_model.lock
    mkfifo model_to_filter.lock filter_to_model.lock

    set filterhome = ~/.filter$$
    if ( ! -e $filterhome) mkdir $filterhome

    # this starts filter but also returns control back to
    # this script immediately.

    (setenv HOME $filterhome; ${MPI} ./filter) &

    while ( -e filter_to_model.lock )

        set todo=`cat < filter_to_model.lock`
        echo "todo received, value = ${todo}"

        if ( "${todo}" == "finished" ) then
          echo "main script: filter done."
          wait
          break

        else if ( "${todo}" == "advance" ) then

          # the second number below must match the number
          # of ensembles. Also, in input.nml, the advance model
          # command must have -np N with N equal to the number
          # of processors this job is using.

          echo "calling model advance now:"
          ./advance_model.csh 0 ${ensemble_size} filter_control00000 || exit 9

          echo "restarting filter."
          ${MPI} ./wakeup_filter

        else

          echo "main script: unexpected value received."
          break

        endif

    end

    echo "filter finished, removing pipes."
    \rm -f model_to_filter.lock filter_to_model.lock

    if ( -d $filterhome) rmdir $filterhome

endif

#-----------------------------------------------------------------------------
# Move the output to storage after filter completes.
# At this point, all the restart,diagnostic files are in the CENTRALDIR
# and need to be moved to the 'experiment permanent' directory.
# We have had problems with some, but not all, files being moved
# correctly, so we are adding bulletproofing to check to ensure the filesystem
# has completed writing the files, etc. Sometimes we get here before
# all the files have finished being written.
#-----------------------------------------------------------------------------

echo "${JOBNAME} ($JOBID) finished at "`date`
echo "Listing contents of CENTRALDIR before archiving"
ls -l

exit 0

${MOVE} *.data *.meta         ${experiment}/POP
${MOVE} data data.cal         ${experiment}/POP
${MOVE} STD*                  ${experiment}/POP

${MOVE} filter_restart*            ${experiment}/DART
${MOVE} assim_model_state_ud[1-9]* ${experiment}/DART
${MOVE} assim_model_state_ic[1-9]* ${experiment}/DART
${MOVE} analysis.nc                ${experiment}/DART
${MOVE} preassim.nc                ${experiment}/DART
${MOVE} obs_seq.final              ${experiment}/DART
${MOVE} dart_log.out               ${experiment}/DART

# Good style dictates that you save the scripts so you can see what worked.

${COPY} input.nml                  ${experiment}/DART
${COPY} *.csh                      ${experiment}/DART
${COPY} $myname                    ${experiment}/DART

ls -lrt

exit 0


