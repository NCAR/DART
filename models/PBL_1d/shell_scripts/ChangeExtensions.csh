#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

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
else

   # need to check on the extension ... it can only be F or f90 ish ...

   foreach FILE ( module_* driver* *_mod.* tridiag.* )  # loops through all the files

      # rename the file

      set BASE = $FILE:r
      set NEWNAME = ${BASE}.F90

      echo -n "Renaming $FILE to $NEWNAME ... "

      mv $FILE ${NEWNAME} || exit 1

      # now change the path_names files with SED

      foreach PATHNAMEFILE ( ../work/path_names* ) 
         set STRING = "1,$ s#src/$FILE#src/$NEWNAME#"
         sed -e "$STRING" $PATHNAMEFILE >! foo
         mv foo $PATHNAMEFILE
      end

      echo "done."

   end

   echo "$SNAME:t complete."

   @ MYSTATUS = 0  
endif
