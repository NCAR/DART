#!/bin/bash -l

#SBATCH --job-name=test
#SBATCH --account=k1418
#SBATCH -N 4
#SBATCH --error=test.%J.err
#SBATCH --output=test.%J.out
#SBATCH --time=00:10:00
#SBATCH --partition=debug
module swap PrgEnv-cray PrgEnv-intel
cp ../../inputs/2obs_obs_seq.txt obs_seq.txt
cp ../../inputs/input.nml .
#cp full_obs_seq.txt obs_seq.txt
time ../create_ocean_obs
rm -f input.nml obs_seq.txt
