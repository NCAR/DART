#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source: /home/thoar/CVS.REPOS/DART/models/cam/shell_scripts/check_model.csh,v $
# $Name:  $

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
