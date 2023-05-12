#!/usr/bin/env bash

comp_name=OCN

# --- info needed from CESM -----
# This info is avaiable in the python case object
#  it would be good to have the case passed to _do_assimilate
#  in case_run.py

exeroot=$(./xmlquery EXEROOT --value) 
rundir=$(./xmlquery RUNDIR --value)
assimilate=$(./xmlquery DATA_ASSIMILATION_$comp_name --value)
ensemble_size=$(./xmlquery NINST_$comp_name --value)

#-------------------------------

# set filter input.nml options

# run filter
if [ $assimilate = TRUE ]; then
   echo "hello I am running filter"
   #mpirun $(rundir)/filter
fi

