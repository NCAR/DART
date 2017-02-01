#!/bin/ksh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

#-- This is just a wrapper to launch num_states instances
#-- of model advance, here filter is serial but launches
#-- num_states number of model advancement process concurrently

process=$1
num_states=$2
control_file=$3

rm -f assim_model_state_ud.*

if [ ${process} == "0" ] ; then
   source ./env.lsf
   mpirun.lsf ./adv_model_par.csh ${num_states} ${control_file}
fi

\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

