#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN CAM_NO_ASSIMILATE"

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   LINK = '/usr/local/bin/ln -fvs'
   breaksw

   default:
      set   LINK = 'ln -fvs'
   breaksw
endsw

# In CESM1_4 xmlquery must be executed in $CASEROOT.
# This script is executed in in $CASEROOT, so we're set.
setenv CASEROOT       `./xmlquery CASEROOT -value`
setenv CASE           $CASEROOT:t

cd ${CASEROOT}
setenv ensemble_size  `./xmlquery NINST_ATM -value`
setenv archive        `./xmlquery DOUT_S_ROOT -value`
setenv SAVE_NTH_SETS  `./xmlquery DOUT_S_SAVE_EVERY_NTH_RESTART_FILE_SET -value`
cd -

#-------------------------------------------------------------------------
# Preliminary clean up, which can run in the background.
# ATM_forcXX, CESM1_5's new archiver has a mechanism for removing restart file sets 
# which we don't need, but it runs ony after the (multicycle) job finishes.  
# It moves the files which need to be archived from RUNDIR to ../../${CASE}.locked 
# and makes (redundant) HARD LINKS in archive/rest/$date and archive/${component}/rest.
# Those hard links are different names for the same physical file.
# Removing one of them, or the original physical file name, does not make the file 
# disappear; the file can still be accessed by the remaining names.
# So we need to remove all names of a file.
# This will remove them during a multicycle run, reducing the need to stop
# to clear scratch space, and preventing {st,lt}_archive from archiving them.
#-------------------------------------------------------------------------
set archive_locked = $archive:h:h/${CASE}.locked/archive

# Remove some log files
$REMOVE        ${archive}/{atm,dart,ice,lnd,rof}/logs/*_00[0-9][02-9].log.* &
$REMOVE ${archive_locked}/{atm,dart,ice,lnd,rof}/logs/*_00[0-9][02-9].log.* &

# We don't need any of the hard links in the $component/rest directories.
# This assumes that cam*.i.* files are archived into the restart directory
# rather than the init directory.  Check that this is set correctly in env_archive.xml.
$REMOVE        ${archive}/{atm,cpl,ice,lnd,rof}/rest/* &
$REMOVE ${archive_locked}/{atm,cpl,ice,lnd,rof}/rest/* &

# Selectively remove restart sets stored in 'top level' restart directory (LINKED and LOCKED).
@ save_mth_days = $SAVE_NTH_SETS / 4
set dates = `ls -dt ${archive}/rest/*`
# Keep the latest one by starting at 2
set d = 2
while ($d <= $#dates)
   set date = $dates[$d]:t 
   set date_parts = `echo $date | sed -e "s/-/ /g"`
   if ($#date_parts == 4) then
      set day_o_month = $date_parts[3]
      set hour_o_day  = $date_parts[4]
      if ( $hour_o_day !~ "00000"  || \
          ($hour_o_day ~= '00000' && $day_o_month%$save_mth_days != 0) ) then
         $REMOVE        ${archive}/$date &
         $REMOVE ${archive_locked}/$date &
      endif
   endif
   @ d++
end

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.cam_${ensemble_member}.i.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

if ( $ensemble_size == 1 ) then
   set FILE = `head -n 1 rpointer.atm`
else
   set FILE = `head -n 1 rpointer.atm_0001`
endif

set FILE = $FILE:r
set ATM_DATE_EXT = `echo $FILE:e`

#=========================================================================
# As implemented, the input filenames are static in the CESM namelists.
# We must link the new uniquely-named files to static names so that when
# the short-term archiver 'restores' the CESM files, the links are right.
#=========================================================================

if ( $ensemble_size == 1 ) then
      set inst_string = ''
      set ATM_INITIAL_FILENAME = ${CASE}.cam${inst_string}.i.${ATM_DATE_EXT}.nc
      ${LINK} ${ATM_INITIAL_FILENAME}    cam_initial${inst_string}.nc || exit -9
else
   set member = 1
   while ( ${member} <= ${ensemble_size} )
      set inst_string = `printf _%04d $member`
      set ATM_INITIAL_FILENAME = ${CASE}.cam${inst_string}.i.${ATM_DATE_EXT}.nc
      ${LINK} ${ATM_INITIAL_FILENAME}    cam_initial${inst_string}.nc || exit -9
      @ member++
   end
endif

echo "`date` -- END CAM_NO_ASSIMILATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

