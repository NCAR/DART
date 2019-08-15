#!/bin/csh 
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script cleans out the work dir.
#
# If called with no arguments, it assumes the current directory.
# If called with arguments, cd to each one in turn and run the script.
#

set origdir = `pwd`

if ( $#argv < 1 ) then
   set arg = ( . )
else
   set arg = ( $argv[*] )
endif

while ( $#arg >= 1 ) 
  set targetdir = $arg[1]

  shift arg
  cd $targetdir
  echo cleaning $targetdir

  @ has_mkmf = `ls mkmf_* | wc -l` >& /dev/null
  
  if ( $has_mkmf > 0 ) then
    foreach TARGET ( mkmf_* )
  
      csh $TARGET || exit $?
      make clean  || exit $?
  
    end
    \rm Makefile 
    \rm -f .cppdefs
  endif
  
  @ has_cdl = `ls *.cdl | wc -l` >& /dev/null
  
  if ( $has_cdl > 0 ) then
    foreach TARGET ( *.cdl )
  
      set OUTNAME = `basename $TARGET .cdl`.nc
      rm $OUTNAME
  
    end
  endif
  
  \rm -f *.o *.mod 
  \rm -f input.nml.*_default
  \rm -f dart_log.*
  
  cd $origdir
end

exit 0

# <next few lines under version control, do not edit>
# $URL: https://proxy.subversion.ucar.edu/DAReS/DART/trunk/models/wrf/work/quickbuild.csh $
# $Revision: 6256 $
# $Date: 2013-06-12 10:19:10 -0600 (Wed, 12 Jun 2013) $

