#!/bin/bash -l

#SBATCH --job-name=compile_all
#SBATCH --account=k1418
#SBATCH --ntasks=1
#SBATCH --error=compile_all.%J.err
#SBATCH --output=compile_all.%J.out
#SBATCH --time=00:15:00
#SBATCH --partition=debug

cp ../../../build_templates/nodes_executables_mkmf.template ../../../build_templates/mkmf.template
time ./quickbuild.csh -mpi

cp ../../../build_templates/front_node_executables_mkmf.template ../../../build_templates/mkmf.template
time ./quickbuild_front_node_executables.csh

