#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   add_path_name.sh
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Adds a file name to all entries with path_names at the beginning -
# use this for adding new files into the DART mkmf framework.

for pathfile in `ls -1 path_names*`
do
  echo "$1" >> $pathfile
done

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

