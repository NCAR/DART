#!/bin/bash -l

#SBATCH --job-name=test
#SBATCH --account=k1418
#SBATCH -N 4
#SBATCH --error=test.%J.err
#SBATCH --output=test.%J.out
#SBATCH --time=00:10:00
#SBATCH --partition=debug
export FILENV=.assign
assign -R
assign -N swap_endian g:all
cp ../../inputs/input.nml ../../inputs/data ../../inputs/data.cal .
ls filter_ics.000[1-4] > input_list.txt
sed 's/_ics/_restart/g' input_list.txt > output_list.txt
time srun --hint=nomultithread --ntasks-per-node=32 --ntasks-per-socket=16 --cpus-per-task=1 --ntasks-per-core=1 --mem-bind=v,local --cpu-bind=threads ../filter	|| exit 3
