#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This is intended to be run as part of the 'test_dart.csh' framework.
#
# This is named 'workshop_setup.csh' to exploit the logic in
# models/run_tests.csh so that these commands are run in this order.
# The CLM files and obs_seq.in are part of the test tarfile.
# models/run_tests.csh replaces the input.nml with input.nml.testing
#
# TODO  inherit the MPI command from run_tests.csh
#----------------------------------------------------------------------

set MODEL = "CLM"

set MPICMD='mpiexec_mpt'

csh quickbuild.csh

cp clm_dart.clm2.r.2013-07-02-00000.nc  clm_restart.nc
cp clm_dart.clm2.h0.2013-07-02-00000.nc clm_history.nc
cp clm_dart.clm2.h2.2013-07-02-00000.nc clm_vector_history.nc

echo 'running clm_to_dart'
./clm_to_dart                    || exit 1

echo 'running perfect_model_obs'
${MPICMD} ./perfect_model_obs    || exit 2

echo 'running fill_inflation_restart'
./fill_inflation_restart         || exit 3

# We need an ensemble of 5 for this test
# We are perturbing a single instance, which will issue a warning.
# At this point, clm_restart.nc has had all the indeterminate values
# replaced and is being copied to all the output files.
# The ensemble is being created by perturbing the single input state.
# We know that will issue a warning, and that's OK.

@ member = 1
while ($member <= 5)
   set F1 = `printf clm_restart_%04d.nc $member`
   set F2 = `printf clm_history_%04d.nc $member`
   set F3 = `printf clm_vector_history_%04d.nc  $member`
   cp clm_restart.nc        $F1
   cp clm_history.nc        $F2
   cp clm_vector_history.nc $F3
   @ member ++
end

ls -1 clm_restart_????.nc        >! restart_files.txt
ls -1 clm_history_????.nc        >! history_files.txt
ls -1 clm_vector_history_????.nc >! vector_files.txt

echo 'running filter'
${MPICMD} ./filter               || exit 4

# For testing purposes, we only need to run dart_to_clm on one
# instance.

echo 'running dart_to_clm'
./dart_to_clm                    || exit 5

exit 0

