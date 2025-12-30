#!/bin/bash

paramfile="$1"

# Load parameter file
source "$paramfile"

# Move to ICBC directory
cd "$ICBC_DIR" || { echo "ERROR: cannot cd to $ICBC_DIR"; exit 1; }

# Run real.exe with MPI
#mpiexec -n 128 -ppn 128 "${RUN_DIR}/WRF_RUN/real.exe"
mpiexec -n 4 -ppn 4 "${RUN_DIR}/WRF_RUN/real.exe"

# Unconditionally set real_done, matching the last two lines of original script
touch "${ICBC_DIR}/real_done"

exit 0

