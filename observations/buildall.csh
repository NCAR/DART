#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section, 
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

set SNAME = $0
set clobber


foreach dir ( */work )
 cd $dir
 echo building in $dir
 echo

 ./quickbuild.csh || exit 1

 cd ../..
end

exit 0
