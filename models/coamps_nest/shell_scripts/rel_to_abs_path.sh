#!/bin/bash
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Convert a file's relative location to an absolute location.
#
# This is taken from
# http://www.shelldorado.com/shelltips/script_programmer.html
#
# The command line argument is the file's relative location
REL_OUTPUT_FILE=$1
D_NAME=`dirname "$REL_OUTPUT_FILE"`
F_NAME=`basename "$REL_OUTPUT_FILE"`
ABS_PATH=`cd "$D_NAME" 2>/dev/null && pwd || echo \"$D_NAME\"`
OUTPUT_FILE="${ABS_PATH}/${F_NAME}"
echo "$OUTPUT_FILE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

