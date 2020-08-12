#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#-------------------------------------------------------------------------------
# DESCRIPTION:
#
# Since each DART diagnostic file at each timestep has been written to a 
# separate file, it is desirable to consolidate all the timesteps for a 
# single ensemble member into a single file.
#
# This must be run AFTER the 'time' dimension has been added to the variables.
# DART_cleanup_add_time.csh  does exactly that.
#
# requires 1 argument:
#    $1 - the ensemble member number
#-------------------------------------------------------------------------------

if ($# != 1) then
   echo "usage: $0 <ensemble_member_number>"
   exit -1
endif 

set MEMBER = $1

#===============================================================================
# everyone has time dimension ... just concatenate them.

\rm -f member_${MEMBER}_list.txt

foreach STAGE( input     preassim    analysis     output )

   set MYSTRING = `printf %s_member_%04d $STAGE $MEMBER`

   find output -name "${MYSTRING}*.out.nc" | grep _d01 | sort >! member_${MEMBER}_list.txt
   if ( `cat member_${MEMBER}_list.txt | wc -l` > 0 ) then

      set OUTPUT = ${MYSTRING}_d01.nc
      echo -n "Creating ${OUTPUT} ... "
      cat member_${MEMBER}_list.txt | ncrcat -O -h -H -o ${OUTPUT}
      echo "done."

   else
      echo "No ${MYSTRING}_d01 files to process."
   endif

   find output -name "${MYSTRING}*.out.nc" | grep _d02 | sort >! member_${MEMBER}_list.txt
   if ( `cat member_${MEMBER}_list.txt | wc -l` > 0 ) then

      set OUTPUT = ${MYSTRING}_d02.nc
      echo -n "Creating ${OUTPUT} ... "
      cat member_${MEMBER}_list.txt | ncrcat -O -h -H -o ${OUTPUT}
      echo "done."

   else
      echo "No ${MYSTRING}_d02 files to process."
   endif

   find output -name "${MYSTRING}*.out.nc" | grep -v _d0 | sort >! member_${MEMBER}_list.txt
   if ( `cat member_${MEMBER}_list.txt | wc -l` > 0 ) then

      set OUTPUT = ${MYSTRING}.nc
      echo -n "Creating ${OUTPUT} ... "
      cat member_${MEMBER}_list.txt | ncrcat -O -h -H -o ${OUTPUT}
      echo "done."

   else
      echo "No $MYSTRING files to process."
   endif

   \rm -f `cat member_${MEMBER}_list.txt`
   \rm -f      member_${MEMBER}_list.txt

end

exit 0

#===============================================================================
