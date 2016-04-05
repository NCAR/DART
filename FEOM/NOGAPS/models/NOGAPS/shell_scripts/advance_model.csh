#!/bin/csh
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
# 1) stages the files into the directory in which to run a model instance
# 2) copies/converts the DART state vector to something the model can ingest
# 3) runs the model
# 4) copies/converts the model output to input expected by DART

set      process = $1
set   num_states = $2
set control_file = $3

# record the DART filter directory (i.e. CENTRALDIR), before cd'ing 
# into the individual directory for the nogaps working directories.
# CENTRALDIR is where the ic and ud files will be,
# as well as where the filter_controlxxxxx file is.

set CENTRALDIR = `pwd`
echo "advance_model.csh running on host "`hostname`

# source an ascii script that sets a bunch of env vars for things
# like top dir, scratch dir, nogaps exec dir, etc.
# FIXME: scratch_dir sounds too temporary for something that's really
# a working dir and needs to persist between invocations of filter.

source ./config.csh

#-------------------------------------------------------------------------
# loop over each state/ensemble member in the filter_control file.
# The filenames in the control file are relative to CENTRALDIR.
#-------------------------------------------------------------------------

set           state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)

   set ensemble_member = `head -$ensemble_member_line $control_file | tail -1`
   set input_file      = `head -$input_file_line      $control_file | tail -1`
   set output_file     = `head -$output_file_line     $control_file | tail -1`

   # Get unique name for the working directory for this member
   set dirname = `printf "mbr%03d" $ensemble_member`

   #----------------------------------------------------------------------
   # Block 0: Error Checking
   #          The ensemble member directories should exist and be ready 
   #          to use because they were made by run_nogapsIC_to_dart.csh 
   #----------------------------------------------------------------------

   set dirstat = 0
   if ( ! -e ${dirname}/bin ) then
      echo "ERROR : ${dirname}/bin does not exist in "`pwd`
      dirstat = 1
   endif

   if ( ! -e ${dirname}/log ) then
      echo "ERROR : ${dirname}/log does not exist in "`pwd`
      set dirstat = 1
   endif

   if ( ! -e ${dirname}/outp ) then
      echo "ERROR : ${dirname}/outp does not exist in "`pwd`
      set dirstat = 1
   endif

   if ( ! -e ${dirname}/specfiles ) then
      echo "ERROR : ${dirname}/specfiles does not exist in "`pwd`
      set dirstat = 1
   endif

   if ( dirstat == 1 ) then
      echo "ERROR : run_nogapsIC_to_dart.csh was supposed to make one."
      exit 1
   endif

   #----------------------------------------------------------------------
   # Block 1: populate a run-time directory with the required bits 
   #----------------------------------------------------------------------

   cd ${dirname}

   # The executables we will need can safely be assumed to exist
   # in the experiment directory (AKA CENTRALDIR) (i.e. ../)

   cp ../input.nml .

   # Run trans_time on the input file to get the target time info.
   # trans_time reads the time from the start of the DART restart file
   # and extracts this time info *without* the overhead of going through
   # dart_to_nogaps -- which requires running part of NOGAPS -- which
   # requires knowledge of 'dtg' -- which is what we are trying to find!
   # i.e. dart_to_nogaps requires a CRDATE.dat file; trans_time does not.

   # trans_time creates "time_info" with 3 lines: now, advance, and delta hours
   ln -sf ../$input_file  temp_ic
   ../${TRANSTIME_exec_name:t} || exit 1
   \rm -f temp_ic

   set     dtg = `head -1 time_info`
   set dtgnext = `head -2 time_info | tail -1`
   set  endtau = `head -3 time_info | tail -1`
   echo "ensemble member $ensemble_member thinks valid date  is $dtg"
   echo "ensemble member $ensemble_member thinks future date is $dtgnext"
   echo "ensemble member $ensemble_member thinks endtau      is $endtau"

   ##---- make the auxiliary analysis-file archive directory for the next cycle 
   #
   #mkdir -p ${scratch_dir}/${dtgnext}/spect${resolution}/${dirname}

   #---- Copy auxiliary analysis files for $dtg (e.g. c3grid, c3land, trpfil, noggeom.txt)
   #---- NOTE:  The shist000000 file copied below is an impostor, used only to read dsqgeot
   #----        The actual shist000000 file is copied after dart->nogaps 

   #---- FIXME: at some point we need to rename/move the $taue version 
   #----        into 000000 - maybe at the end of this run?
   #---- TJH asks ... are these copies backwards?
   # cp ${scratch_dir}/${dtg}/spect${resolution}/${dirname}/shist000000  \
   #                                   ${dirname}/specfiles/shist000000 
   # cp ${scratch_dir}/${dtg}/spect${resolution}/${dirname}/c3grid000000 \
   #                                   ${dirname}/specfiles/c3grid000000
   # cp ${scratch_dir}/${dtg}/spect${resolution}/${dirname}/c3land000000 \
   #                                   ${dirname}/specfiles/c3land000000
   # cp ${scratch_dir}/${dtg}/spect${resolution}/${dirname}/trpfil       \
   #                                   ${dirname}/specfiles/trpfil
   # cp ${scratch_dir}/climo${resolution}/noggeom.txt                    \
   #                 ${dirname}/specfiles/noggeom.txt


   # we assume that taui (initial time) is always going to be 0
   # (the restart file is at the current time).   we believe this,
   # and having only one time level in the restart file, will make
   # the model do a euler forward time step before starting the two
   # level leapfrog scheme.  

   # if taui is 0, the restart files always have to be named shist000000
   # so either dart_to_nogaps should just output this as the default
   # output name, or this script will need to rename that file here.

   # might need to update tauh - how often to output history file if
   # advancing more than a single time period.  can set = to taue?

cat <<EOF3 > rdifil
 &rdilst
 iproc  = $iproc,
 jproc  = $jproc,
 jsplit = $jsplit,
 lgtrdy = f,
 lnmode = t,
 loutp  = f,
 lfcst  = t,
 lmpi2  = f,
 al2al  = t/
EOF3

cat <<EOF4 > filist
 &namfil
 ocards  = 'isisocd',
 idfile  = 'idisis',
 ifilin  = './specfiles',
 ifilout = './outp',
 sstdir  = './f/f$1',
 icedir  = './f/f$1',
 snodir  = './f/f$1',
 pstdir0 = './',
 pstdir1 = './specfiles',
 hstdir  = './specfiles',
 clmdir  = '${climo}'/
EOF4

cat <<EOF5 > namlsts
 &modlst
 dt      = 200.,
 dtrad   =  12.0,
 taui    =   0.0,
 taue    = $endtau,
 tauh    = $endtau,
 vistsh  = .0000354,.000707,2.828,
 visvor  = .0000354,.000707,2.828,
 visdiv  = .0000354,.000707,2.828,
 dtgfnoc = '${dtg}',
 nmodev  = 3,
 jskip   = 2/
 &keylst
 lgeosfilt=t/
EOF5

   cp ${ocards_files}/idisis.th.car                idisis
   cp ${ocards_files}/isisocd.jun05jul05.32.scaled isisocd

   #----------------------------------------------------------------------
   # Block 2: Convert DART file to form needed by NOGAPS      dart->nogaps
   #----------------------------------------------------------------------
   # dart_to_nogaps expects the $input_file to be a certain name.
   # input.nml:&dart_to_nogaps_nml:dart_to_nogaps_input_file = 'dart_old_vector',
   #
   # dart_to_nogaps creates a three-line file containing the header info from 
   # the dart initial conditions file (namelist default is 'dart_data.time')
   # which is completely superfluous as far as advance_model.csh 
   # is concerned.

   mv ../$input_file        dart_old_vector

   # CRDATE.dat is the current time - read by NOGAPS codes
   echo     $dtg >! CRDATE.dat

   ${MPI} ../${DART2NOGAPS_exec_name:t}

   if ( $status != 0 ) then
          echo "status = $status"
          echo "ERROR: dart to nogaps conversion failed for member $ensemble_member"
          echo "ERROR: check contents of "`pwd`
          exit 2
   endif

   #---- position the $dtg analysis shist file for archival 
   #
   # FIXME: if you want to save the restarts for later forecasts,
   #  copy this shist000000 to something with dtg in it.
   #cp  ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist000000 \
   #    ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist.${dtg}

   #cp  ${scratch_dir}/${experiment_name}/${dirname}/specfiles/shist000000 \
   #    ${scratch_dir}/${dtg}/spect${resolution}/${dirname}

   #----------------------------------------------------------------------
   # Block 3: advance the model
   #----------------------------------------------------------------------

   ${MPI} ../${NOGAPS_exec_name}

   if ( $status != 0 ) then
          echo "ERROR: status = $status"
          echo "ERROR: nogaps forecast failed for member $ensemble_member at "`cat CRDATE.dat`
          exit 3
   endif

   #---- copy the taue forecasts to the $dtgnext auxiliary analysis-file archive directory

   set leadtime  = `printf %06d $endtau`

   # FIXME: if you want to save the output from the nogaps run, use
   # something like:
   #cp specfiles/shist${leadtime} ${someplace}/${dirname}/shist.${dtgnext}

   mv specfiles/shist$leadtime  specfiles/shist000000
   mv specfiles/c3grid$leadtime specfiles/c3grid000000
   mv specfiles/c3land$leadtime specfiles/c3land000000

   #----------------------------------------------------------------------
   # Block 4: Convert NOGAPS state to a DART vector: nogaps->dart
   #----------------------------------------------------------------------
   # &nogaps_to_dart_nml:nogaps_to_dart_output_file must have
   # the value 'dart_new_vector' for this logic to work.
   #---------------------------------------------------------

   echo $dtgnext >! CRDATE.dat
   
   ${MPI} ../${NOGAPS2DART_exec_name:t}

   if ( $status != 0 ) then
          echo "ERROR: status = $status"
          echo "ERROR: nogaps to dart conversion failed for member $ensemble_member at "`cat CRDATE.dat`
          exit 4
   endif

   #---- move the new dart state vector to where filter 
   #---- expects it to be.

   mv dart_new_vector ../$output_file || exit 4

   #---- move back in position for the next ensemble member

   cd ..

   #---- increment counters for next loop

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line      = $input_file_line      + 3
   @ output_file_line     = $output_file_line     + 3

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


#---- MANDATORY - Remove the control_file to signal completion. If it still
#---- exists in CENTRALDIR after all the ensemble members have been advanced,
#---- it means one or more of the advances failed and is an ERROR CONDITION.

\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

# original collaboration information - do not edit
# $orgURL: https://svn2.assembla.com/svn/ngdart/shell_scripts/advance_model.csh $
# $orgId: advance_model.csh 111 2010-06-09 21:55:44Z thoar $
# $orgRevision: 111 $
# $orgDate: 2010-06-09 15:55:44 -0600 (Wed, 09 Jun 2010) $

