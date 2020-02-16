#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Display the destination files from repack_st_arch.csh.
# This can't be rolled into repack_st_arch.csh
# because that submits globus transfer jobs
# that take an unknown amount of time to finish.

# On the data-access node.

if ($#argv == 0) then
   set y = 2017
   set m = 1
else if ($#argv == 2) then
   set y = $1
   set m = $2
else
   echo "Usage: after running case.st_archive and repack_st_archive,"
   echo "       on casper:"
   echo "       % pre_purge.csh yyyy m"
   exit
endif
set ym = `printf %s-%02d $y $m`

cd $pr/$casename

echo "Coupler history (forcing) files:"
ls -lt cpl/hist/0080/*${y}*

echo "\n Component history files:"
ls -lt {lnd,atm,ice,rof}/hist/0080/*.{clm2,cam,cice,mosart}_*.h*${y}*[cz]

echo "\n DART obs space diagnostic files:"
ls -lt esp/hist/${ym}

echo "\n Restart files in ${csdir}/${casename}/rest/${ym}:"
ls -l ${csdir}/${casename}/rest/${ym}/*{0080,inf}*

echo "\n Ensemble means and inflation files in ${csdir}/${casename}/esp/hist/${ym}:"
ls -l ${csdir}/${casename}/esp/hist/${ym}

echo "\n CAM preassim ensembles in ${csdir}/${casename}/atm/hist/${ym}:"
ls -l ${csdir}/${casename}/atm/hist/${ym}

echo "\n DART log files in ${csdir}/${casename}/logs/${ym}:"
ls -l ${csdir}/${casename}/logs/${ym}

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
