#!/bin/csh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This script cleans up a DART work dir.  It removes all the executables
# built by the quickbuild.csh script, plus other files created by the
# quickbuild.csh script (e.g. Makefile, input.nml.*_default, etc),
# and any .nc files created from .cdl files.  It does *not* remove
# files created by running the executables (e.g. output .nc files,
# obs_seq.* files, etc) other than the dart_log.* files.
#
#
# Usage: 
#  If called with no arguments, it runs in the current directory.
#  If called with arguments, cd to each one in turn and run the script.
#
# Example usage to clean up entire dart checked out tree:
#
#  - put quickclean.csh somewhere on your search path or
#    add $DART/assimilation_code/scripts to your path.
#
#   > cd $DART
#   > find . -name work -exec quickclean.csh {} \;
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
# $URL:$
# $Revision:$
# $Date:$

