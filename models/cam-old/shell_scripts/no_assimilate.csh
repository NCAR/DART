#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

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

set ensemble_size = ${NINST_ATM}

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


