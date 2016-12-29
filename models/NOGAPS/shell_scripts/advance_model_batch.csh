#!/bin/csh -v
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# NOGAPS VERSION
#
# Arguments are (created by 'filter' or 'perfect_model_obs' and include):
# 1) the process number of caller,
# 2) the number of ensemble members/state copies belonging to that process, and 
# 3) the name of the control_file for that process.
# 
# If this script finishes and the 'control_file' still exists, it is
# an ERROR CONDITION and means one or more of the ensemble members did
# not advance properly. Despite our best attempts to trap on this
# condition, some MPI installations simply hang, some properly terminate.
#
# This script loops over all the entries in the control_file to advance 
# any/all of the ensemble members.  The number of trips through the 
# loop is the second argument to this script. The control_file contains 
# the information about which ensemble members are to be advanced by THIS TASK.
# Sometimes it may be just one ensemble member, sometimes all of them.
# Read DART/doc/html/filter_async_modes.html and the mpi_intro.html
# for an overview.
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance
#    and copies the necessary files into the temporary directory
# 2) copies/converts the DART state vector to something the model can ingest
# 3) runs the model
# 4) copies/converts the model output to input expected by DART

# this is the batch version - it won't be called with arguments.
# these need values.  the first is always 0, the last is always
# filter_control000000 - the second one must be the number of
# ensembles.  for now, put them on the command line when you
# call this script.

set      process = $1
set   num_states = $2
set control_file = $3

echo 'in advance model'

# FIXME: scratch_dir sounds too temporary for something that's really
# a working dir and needs to persist between invocations of filter.

# source an ascii .csh script that sets a bunch of env vars for things
# like top dir, scratch dir, nogaps exec dir, etc.

source ./config.csh

echo "changing directories to ${scratch_dir}/${experiment_name}"
cd ${scratch_dir}/${experiment_name}

# update the filter_restart files from the previous step to
# be model_advance files (with 2 timestamps - the go-to time
# and the current time.

# FIXME: we need ADVDAY and ADVSEC from $dtgnext?
# something like set ADV=(`echo $dtgnext +0 -g | ./advance_time`)
# create restart.nml with ADVDAY and ADVSEC and just the &restart_file_tool
# name list in it - and remove &restart_file_tool from default input.nml
# set ADVDAY = $ADV[1]
# set ADVSEC = $ADV[2]
# sed -e "s/ADVDAY/$ADVAY/" \
#     -e "s/ADVSEC/$ADVSEC/" < restart.nml >> input.nml
#./${RESTART_exec_name:t}

# Determine the number of ensemble members from input.nml,

# pull out the first input_file line from the control file here,
# outside the loop, run trans_time just to get the time info.
# it reads the time from the start of the restart/model_advance 
# file and makes a time file without doing anything to the restart data.
# endtau we will know after running trans_time. 

set input_file = `head -2 $control_file | tail -1`
echo $input_file
pwd

# copy this file to a known filename so we can run trans_time.
# this creates time_info with 3 lines: now, advance, and delta hours
ln -s $input_file  temp_ic
./trans_time
rm -f temp_ic

# the first line is the dtg time now, which is what we need
set time_now = `head -1 ./time_info`
set dtg = $time_now
echo date is $dtg

set endtau = `head -3 ./time_info | tail -1`

# at the end, we can use this time to store the updated
# data - use this dtgnext time to name the directories for
# where the files are going to go.
# FIXME: this executable should be in assembla and built when
# the executables are built.
set dtgnext = `echo $dtg\ +$endtau | /home/coral/hansenj/bin/newdtg`
echo dtgnext is $dtgnext


##---- make the auxiliary analysis-file archive directory for the next cycle 
#
#mkdir -p ${scratch_dir}/$dtgnext
#mkdir -p ${scratch_dir}/${dtgnext}/spect$resolution
#

#---- make the assim_model_state_ud.* archive directory for the next cycle

# ??
#mkdir -p ${scratch_dir}/${dtgnext}/ud



#########################################################
#
#  loop through each state/ensemble member in the
#  filter_control file.
#
#########################################################


set state_copy           = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line $control_file | tail -1`
   set input_file      = `head -$input_file_line      $control_file | tail -1`
   set output_file     = `head -$output_file_line     $control_file | tail -1`
   
   set index   = `printf %03d $ensemble_member`
   set dirname = 'mbr'${index}

#---- the ensemble member directories in the scratch directory should exist
#---- and be ready to use.

   #mkdir -p $dirname/bin
   #mkdir -p $dirname/log
   #mkdir -p $dirname/specfiles
   #mkdir -p $dirname/outp

#---- make the given member's auxiliary analysis-file archive directory for the next cycle

#   mkdir -p ${scratch_dir}/${dtgnext}/spect${resolution}/$dirname

#---- Copy auxiliary analysis files for $dtg (e.g. c3grid, c3land, trpfil, noggeom.txt)
#---- NOTE:  The shist000000 file copied below is an impostor, used only to read dsqgeot
#----        The actual shist000000 file is copied after dart->nogaps 

# FIXME: at some point we need to rename/move the $taue version of these
# into 000000 - maybe at the end of this run?
#   cp ${scratch_dir}/$dtg/spect${resolution}/$dirname/shist000000              $dirname/specfiles/shist000000 
#   cp ${scratch_dir}/$dtg/spect${resolution}/$dirname/c3grid000000             $dirname/specfiles/c3grid000000
#   cp ${scratch_dir}/$dtg/spect${resolution}/$dirname/c3land000000             $dirname/specfiles/c3land000000
#   cp ${scratch_dir}/$dtg/spect${resolution}/$dirname/trpfil                   $dirname/specfiles/trpfil
#   cp ${scratch_dir}/climo${resolution}/noggeom.txt                            $dirname/specfiles/noggeom.txt

#---- collect the 3 executables we will need in the bin dir.

   cp ${NOGAPS_exec_dir}$NOGAPS_exec_name          ${dirname}/bin/
   cp $DART2NOGAPS_exec_name                       ${dirname}/bin/
   cp $NOGAPS2DART_exec_name                       ${dirname}/bin/
   set NOGAPS_exec_mpi      = ${scratch_dir}/${experiment_name}/${dirname}/bin/$NOGAPS_exec_name
   set DART2NOGAPS_exec_mpi = ${scratch_dir}/${experiment_name}/${dirname}/bin/${DART2NOGAPS_exec_name:t}
   set NOGAPS2DART_exec_mpi = ${scratch_dir}/${experiment_name}/${dirname}/bin/${NOGAPS2DART_exec_name:t}

   # FIXME:  is this right?  (we think yes; it seems to be needed...)
   cp CRDATE.dat                 ${dirname}
   cp input.nml                        ${dirname}

#---- move $input_file to some standard name

   cp $input_file                      ${dirname}/dart_old_vector

# we assume that taui (initial time) is always going to be 0
# (the restart file is at the current time).   we believe this,
# and having only one time level in the restart file, will make
# the model do a euler forward time step before starting the two
# level leapfrog scheme.  

# if taui is 0, the restart files always have to be named shist000000
# so either dart_to_nogaps should just output this as the default
# output name, or this script will need to rename that file here.

#  might need to update tauh - how often to output history file if
#  advancing more than a single time period.  can set = to taue?

   

#################################################    
#
# Now generate a submit script for the given ensemble member
#

   if (-d $dirname) then
        echo member_name is $dirname
   endif

        cat <<EOFINIT > "${dirname}/bin/ufcst$dirname"
#!/bin/csh 
#
#########################################################
#
#    This file is submitted to the LSF schedular
#
#########################################################

# BSUB -J m${index}
# BSUB -q economy
# BSUB -W $clock_time
# BSUB -n $nproc
# BSUB -o ${NOGAPS_log_dir}ufcst$dirname.out
# BSUB -e ${NOGAPS_log_dir}ufcst$dirname.err
# BSUB -K

########################################
#
#  Set some environment variables that will be global to all of the
#  ensemble member runs.
#
########################################

# The home directory
# (location of standard output and standard error files
#  created by the script submitted to the LoadLeveler)
#

  set HOME     = ${scratch_dir}/${experiment_name}/${dirname}

#
# terrain and climotological files
#
  set climo    = ${scratch_dir}/climo$resolution
#
#
# The data directory
# (location of read-only data like initial fields
#

  set DATA=\$HOME/specfiles
  set GFLD=\$HOME/outp

#
# The working scratch directory.
# (where output data is written during the model run)
#

  set TMP=\$HOME

#
# The initial fields directory
#

  set flds=\$HOME/specfiles

#
#
########################################

set echo
echo $dtg >! \$TMP/CRDATE.dat

# create a logical link to \$TMP.  This is sometimes needed
# if you use very long pathnames that exceed the size of the 
# character variables that hold them in the Fortran code.

ln -fs \$TMP temp

######################################## 

cd \$HOME

cat <<EOF3 > rdifil
 &rdilst
 iproc = $iproc,
 jproc = $jproc,
 jsplit= $jsplit,
 lgtrdy=f,
 lnmode=t,
 loutp= f,
 lfcst= t,
 lmpi2=f,
 al2al=t/
EOF3

cat <<EOF4 > filist
 &namfil
 ocards='isisocd',
 idfile='idisis',
 ifilin='\$TMP/specfiles/',
 ifilout='\$TMP/outp/',
 sstdir='temp/f/f$1',
 icedir='temp/f/f$1',
 snodir='temp/f/f$1',
 pstdir0='\$TMP/',
 pstdir1='\$TMP/specfiles/',
 hstdir='\$TMP/specfiles/',
 clmdir='\$climo/'/
EOF4

cat <<EOF5 > namlsts
 &modlst
 dt=200.,
 dtrad= 12.0,
 taui=   0.0,
 taue= $endtau,
 tauh= $endtau,
 vistsh=.0000354,.000707,2.828,
 visvor=.0000354,.000707,2.828,
 visdiv=.0000354,.000707,2.828,
 dtgfnoc="$dtg",
 nmodev=3,
 jskip=2/
 &keylst
 lgeosfilt=t/
EOF5

cp ${ocards_files}/idisis.th.car                idisis
cp ${ocards_files}/isisocd.jun05jul05.32.scaled isisocd



#---------------------------------------------------------
# run dart->nogaps 
#---------------------------------------------------------

$MPI $DART2NOGAPS_exec_mpi

if ( \$status != 0 ) then
       echo status = \$status
       echo "dart to nogaps conversion failed "
       exit
endif

#---- position the $dtg analysis shist file for archival 
#
# FIXME: if you want to save the restarts for later forecasts,
#  copy this shist000000 to something with dtg in it.
#cp  ${scratch_dir}/${experiment_name/${dirname}/specfiles/shist000000 \
#    ${scratch_dir}/${experiment_name/${dirname}/specfiles/shist.${dtg}

#cp  ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist000000 \
#    ${scratch_dir}/${dtg}/spect${resolution}/${dirname}



#-------------------------
# advance the model 
#-------------------------

$MPI $NOGAPS_exec_mpi

if ( \$status != 0 ) then
       echo status = \$status
       echo "nogaps forecast failed "
       exit
endif

#---- copy the taue forecasts to the $dtgnext auxiliary analysis-file archive directory

echo ready to set leadtime here
set leadtime  = `printf %06d $endtau`
echo leadtime now \$leadtime

# FIXME: check that this is doing the right thing.  the original
# code is here, but we want to keep the dtgnext out of the pathname.
#cp ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist\${leadtime} \
#   ${scratch_dir}/${dtgnext}/spect${resolution}/${dirname}/shist000000

# FIXME: if you want to save the output from the nogaps run, use
# something like this.
#cp ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist\$leadtime \
#   ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist.${dtgnext}

mv ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist\${leadtime} \
   ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist000000
mv ${scratch_dir}/${experiment_name}/${dirname}/specfiles/c3grid\${leadtime} \
   ${scratch_dir}/${experiment_name}/${dirname}/specfiles/c3grid000000
mv ${scratch_dir}/${experiment_name}/${dirname}/specfiles/c3land\${leadtime} \
   ${scratch_dir}/${experiment_name}/${dirname}/specfiles/c3land000000

#---------------------------------------------------------
# run nogaps->dart
#---------------------------------------------------------

# echo dtgnext to the file CRDATE.dat (to be used by nogaps_to_dart)
echo $dtgnext >! CRDATE.dat

$MPI $NOGAPS2DART_exec_mpi

if ( \$status != 0 ) then
       echo status = \$status
       echo "nogaps to dart conversion failed "
       exit
endif

#---- move the new dart state vector to where filter is going to expect
#---- to read it.

mv ${scratch_dir}/${experiment_name}/${dirname}/dart_new_vector \
   ${scratch_dir}/${experiment_name}/${output_file}

#endif

EOFINIT



chmod 777 "${dirname}/bin/ufcst$dirname"

if (! -e ${scratch_dir}/${experiment_name}/${dirname}/bin/submit_script) then

    cat <<EOF >> ${scratch_dir}/${experiment_name}/${dirname}/bin/submit_script
#!/bin/csh 
# This script submits all the jobs to the batch system.
#
# Check and see how many jobs are already queued.
set count = `bjobs -u \$USER | fgrep \$USER | wc -l`
echo "There are \$count jobs already in the queue"
EOF
chmod +x ${scratch_dir}/${experiment_name}/${dirname}/bin/submit_script

endif

cat <<EOF >> ${scratch_dir}/${experiment_name}/${dirname}/bin/submit_script

bsub < ${scratch_dir}/${experiment_name}/${dirname}/bin/ufcst$dirname

EOF

#---- Run the batch script to submit the jobs.

   cd ${scratch_dir}/${experiment_name}/${dirname}
   bin/submit_script
   cd ${scratch_dir}/${experiment_name}

#---- move back in position for the next ensemble member

cd ${scratch_dir}/${experiment_name}

#---- increment counters for next loop

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line      = $input_file_line + 3
   @ output_file_line     = $output_file_line + 3


#---------------------------------------
#---- end of loop over ensemble members
#---------------------------------------

end


#---- tar-up the analysis files for $dtg

#cd $scratch_dir
#tar -cf new.anal.${dtg}.tar ${dtg}/spect$resolution
#mv new.anal.${dtg}.tar $archive_analysis_dir 


#---- tar-up the auxiliary analysis files for $dtgnext

#tar -cf aux.anal.${dtgnext}.tar ${dtgnext}/spect$resolution
#mv aux.anal.${dtgnext}.tar $archive_analysis_dir


#---- tar-up the assim_model_state_ud.* files for $dtgnext

#tar -cf assim.model.state.ud.files_${dtgnext}.tar ${dtgnext}/ud
#mv assim.model.state.ud.files_${dtgnext}.tar $archive_analysis_dir 


#---- work-directory clean-up
#---- If you are debugging, you may want to keep these directories. 

#cd $scratch_dir
#\rm -fr $dtg $dtgnext climo$resolution $experiment_name

# FIXME: right now, we are submitting these batch jobs with the -K
# flag which says run this batch job and wait for it to finish.
# we really should submit all of them and then hang out here waiting
# for them to exit or finish.  i wonder if we can start a dummy job
# that only runs when all the ensembles have either run or exited with
# an error, and after it runs, make sure all the output files are done
# otherwise exit here with -1.

#---- MANDATORY - Remove the control_file to signal completion. If it still
#---- exists in CENTRALDIR after all the ensemble members have been advanced,
#---- it means one or more of the advances failed and is an ERROR CONDITION.

#\rm -rf $control_file


exit 0


# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

# original collaboration information - do not edit
# $orgURL: https://svn2.assembla.com/svn/ngdart/shell_scripts/advance_model_batch.csh $
# $orgId: advance_model_batch.csh 110 2010-06-09 21:52:58Z thoar $
# $orgRevision: 110 $
# $orgDate: 2010-06-09 15:52:58 -0600 (Wed, 09 Jun 2010) $
