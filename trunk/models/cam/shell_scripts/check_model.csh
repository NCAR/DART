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

if ($#argv < 1) then
   echo "usage; check_cam num_ens_members"
   exit
endif 

set n = 1
while ($n <= $1)
     tail -40 cam_out_temp$n | grep 'END OF MODEL RUN' > /dev/null
     if ($status != 0) echo cam_out_temp$n finished abnormally
     @ n++
end
exit
