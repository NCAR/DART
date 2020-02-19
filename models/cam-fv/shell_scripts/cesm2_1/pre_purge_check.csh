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
# Usage: after case.st_archive and repack_st_archive are finished,
#        on casper:
#        % pre_purge.csh 

# On the data-access node.

source ./data_scripts.csh

set ym = `printf %s-%02d $data_year $data_month`

cd $data_proj_space

echo "Coupler history (forcing) files:"
ls -lt cpl/hist/00${data_NINST}/*${data_year}*

echo "\n Component history files:"
ls -lt {lnd,atm,ice,rof}/hist/00${data_NINST}/*.{clm2,cam,cice,mosart}_*.h*${data_year}*[cz]

echo "\n DART obs space diagnostic files:"
ls -lt esp/hist/${ym}

echo "\n Restart files in ${data_campaign}/${data_CASE}/rest/${ym}:"
ls -l ${data_campaign}/${data_CASE}/rest/${ym}/*{00${data_NINST},inf}*

echo "\n Ensemble means and inflation files in ${data_campaign}/${data_CASE}/esp/hist/${ym}:"
ls -l ${data_campaign}/${data_CASE}/esp/hist/${ym}

echo "\n CAM preassim ensembles in ${data_campaign}/${data_CASE}/atm/hist/${ym}:"
ls -l ${data_campaign}/${data_CASE}/atm/hist/${ym}

echo "\n DART log files in ${data_campaign}/${data_CASE}/logs/${ym}:"
ls -l ${data_campaign}/${data_CASE}/logs/${ym}

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
