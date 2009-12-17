#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This script renames the PBL_1d source files from the (nonstandard)
# .F (implying f77 syntax and preprocessing) to .f90 (reflecting the actual syntax)

set noclobber
set SNAME = $0

if ( $cwd:t != 'src' ) then
   echo
   echo "ERROR: Must execute this command from the PBL_1d/src directory."
   echo "(need to know relative position so we can change path_names files)"
   echo
   @ MYSTATUS = 2
   exit
endif

# Echo usage notes _or_ start
 
if ($#argv != 1) then
   echo "usage: $SNAME:t DesiredExtension"
   echo "This script changes the extension on every file in PBL_1d/src"
   echo "and modifies the path_names_xxxx files accordingly."
   echo
   echo "for example, if you want all the files to end in .f90"
   echo "the appropriate command is "
   echo "$SNAME:t f90"
   echo
   @ MYSTATUS = 1  
   exit
endif

# check to see if the requested extension is reasonable

switch ( $1 )
case F90:
      set FEXT = $1
      breaksw
case F:
      set FEXT = $1
      breaksw
case f90:
      set FEXT = $1
      breaksw
case f:
      set FEXT = $1
      breaksw
default:
      echo "unrecognized extension ( $1 ) doing nothing"
      exit
endsw

# loop through all the files

foreach FILE ( module_* driver* *_mod.* tridiag.* )

   # rename the file

   set BASE = $FILE:r
   set NEWNAME = ${BASE}.${FEXT}

   echo "Renaming $FILE to $NEWNAME ... "

   mv $FILE ${NEWNAME} || exit 1
   # svn rename $FILE ${NEWNAME} || exit 1

   # now change the path_names files with SED

   foreach PATHNAMEFILE ( ../work/path_names* ) 
      set STRING = "1,$ s#src/$FILE#src/$NEWNAME#"
      sed -e "$STRING" $PATHNAMEFILE >! foo
      mv foo $PATHNAMEFILE
   end

end

echo "$SNAME:t complete."

@ MYSTATUS = 0  

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

