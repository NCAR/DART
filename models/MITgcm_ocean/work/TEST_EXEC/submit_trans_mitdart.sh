#!/bin/bash -l
#SBATCH --job-name=PROC
#SBATCH --account=k1418
#SBATCH --ntasks=64
#SBATCH --error=PROC.%J.err
#SBATCH --output=PROC.%J.out
#SBATCH --time=30:00
#SBATCH --partition=debug
module swap PrgEnv-cray PrgEnv-intel

rm -rf TEMPORARY
mkdir TEMPORARY
cd TEMPORARY
cp ../../../inputs/input.nml ../../../inputs/data  .
echo "M2D" > MODE.txt	#M2D for MIT2DART, D2M for DART2MIT
for i in $(seq -w 1 0004);do
 rm -f OUTPUT.nc
 ln -f ../../../inputs/MEM_${i}/* .
 ../../trans_mitdart 
 mv OUTPUT.nc ../filter_ics.${i}
done
cd ..
rm -rf TEMPORARY
