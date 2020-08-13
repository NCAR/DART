#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This is an example script for how to stage the files in CENTRALDIR
# in preparation for an assimilation.
#
#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
#==============================================================================

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

#==============================================================================
# Stage all the required files in CENTRALDIR
#
# CENTRALDIR is where 'filter' will run, each model advance takes place
# in a subdirectory created and populated by 'advance_model.csh'
#
# The files may exist from a perfect model setup, use them. If not, get them.
#==============================================================================

set CENTRALDIR = `pwd`

set NOAHDIR = /Users/thoar/svn/DART/devel/models/noah/src/hrldas-v3.3
set DARTDIR = /Users/thoar/svn/DART/devel/models/noah

foreach FILE ( Noah_hrldas_beta SOILPARM.TBL VEGPARM.TBL GENPARM.TBL URBPARM.TBL )
   if ( -e ${FILE} )  then
      echo "Using existing $FILE"
   else
      ${COPY} ${NOAHDIR}/Run/${FILE} . || exit 1
   endif
end

foreach FILE ( wrfinput namelist.hrldas )
   if ( -e ${FILE} ) then
      echo "Using existing ${FILE}"
   else
      ${COPY} ${DARTDIR}/templates/${FILE}.template ${FILE} || exit 1
   endif
end

foreach FILE ( obs_seq.out input.nml filter dart_to_noah noah_to_dart restart_file_tool )
   if ( -e ${FILE} )  then
      echo "Using existing $FILE"
   else
      echo "$FILE needs to be copied."
      ${COPY} ${DARTDIR}/work/${FILE} . || exit 1
   endif
end

foreach FILE ( run_filter.csh advance_model.csh )
   if ( -e ${FILE} )  then
      echo "Using existing $FILE"
   else
      ${COPY} ${DARTDIR}/shell_scripts/${FILE} . || exit 1
   endif
end

#==============================================================================
# need a set of noah restart files to define the initial ensemble
#
# 1) point to a directory full of noah restart files and pick N of them;
#    The tricky part is that the Time variable in those files is all wrong.
# 2) convert them to DART format;
# 3) make sure the time in each DART file is 'identical'
# 4) the original noah restart files are also needed to start up
#    advance_model.csh for the first model advance.
#
# NOTE : This is for testing the machinery ONLY! If you try to publish a paper
# with this ensemble I will REJECT IT AND EVERY OTHER PAPER IN YOUR CAREER.
#
# COME UP WITH YOUR OWN ENSEMBLE.
#==============================================================================

echo "Make sure the earliest time in the obs_seq.out file is at or after"
echo "the time we are inserting in the initial ensemble."
echo "DART can advance the model states to the observation time."
echo "DART cannot move the model state back in time."

set ENSEMBLESOURCE = /Users/thoar/svn/DART/devel/models/noah/ensemble_source

set   nfiles = `ls -1 ${ENSEMBLESOURCE}/RESTART*DOMAIN* | wc -l` || exit 2
set filelist = `ls -1 ${ENSEMBLESOURCE}/RESTART*DOMAIN*`         || exit 2

# now that we know everything about the initial ensemble; 
# ensure the times are consistent and then convert to DART initial condition files.

@ ifile = 1
@ ensemble_member = 0
while ($ifile <= $nfiles)

   @ ensemble_member = $ensemble_member + 1
   set fext = `printf %04d $ensemble_member`

   # make sure all the initial ensemble files have the same time.
   # If the restart files already have identical, correct times, you can skip this part.

   ncap2 -O -s 'Times(0,:)="2009-01-02_01:00:00"' $filelist[$ifile] restart.$fext.nc
   if ($status != 0) then
      echo "WARNING: time conversion failed"
   endif
   ln -svf restart.$fext.nc restart.nc

   # make initial conditions for DART

   ./noah_to_dart                     || exit 3
   ${MOVE} dart_ics filter_ics.$fext  || exit 4

   @ ifile = $ifile + 1
end

# Since we have some knowledge of the ensemble size, 
# provide reasonable default values.

ex input.nml <<ex_end
/filter
g;ens_size ;s;= .*;= $ensemble_member,;
g;num_output_state_members ;s;= .*;= $ensemble_member,;
g;num_output_obs_members ;s;= .*;= $ensemble_member,;
g;single_restart_file_in ;s;= .*;= .false.,;
g;single_restart_file_out ;s;= .*;= .false.,;
wq
ex_end

# DART needs a copy of the NOAH restart file to determine the sizes
# of the state vector components. We are going to LEAVE the final
# ensemble member restart file linked to the expected restart file name.

#==============================================================================
# Finish up.
#==============================================================================

echo
echo "CENTRALDIR is ${CENTRALDIR}"
echo "Configure     ${CENTRALDIR}/input.nml"
echo "Configure     ${CENTRALDIR}/namelist.hrldas"
echo "Configure     ${CENTRALDIR}/wrfinput"
echo "execute       ${CENTRALDIR}/run_filter.csh"
echo

exit 0


