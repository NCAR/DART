#!/bin/bash

declare -a arr=("gcc" "ifort" "ifx" "cce" "nvhpc")

for compiler in "${arr[@]}"; do

outfile=build_with_$compiler.sh

cat <<EOF > $outfile
#!/bin/bash
#PBS -A P86850054
#PBS -N build-everything-$compiler
#PBS -j oe
#PBS -k eod
#PBS -q main
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=128:mpiprocs=128

export TMPDIR=/tmp/

time ./run_all_quickbuilds.sh $compiler $1

EOF

qsub $outfile

done
