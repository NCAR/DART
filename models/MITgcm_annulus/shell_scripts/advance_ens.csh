#!/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next six lines automatically updated by CVS, do not edit>
# $Source$
# $Revision$
# $Date$
# $Author$
# $Name$ 
# $Id$
#

echo "Entering advance_ens.csh"

setenv PBS_O_WORKDIR `pwd`

### Output to confirm job characteristics
echo "Running on host "`hostname`
echo "Time is "`date`
echo "Directory is "`pwd`

# First line of filter_control should have number of model states to be 
# integrated
set nensmbl = `head -1 filter_control`

set element = 1
while($element <= $nensmbl)
   $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element
   sleep 0.1
   @ element++
end

# Wait for all *background* processes to finish up
wait

# signal to async_filter.csh to continue
rm -f $PBS_O_WORKDIR/batchflag

echo "Leaving advance_ens.csh"
