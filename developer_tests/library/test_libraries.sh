#!/usr/bin/env bash

set -e
export DART=$(git rev-parse --show-toplevel)

main () {

NETCDF_LIB=/opt/local/
testdir=$(pwd)

# remove any previous test directoires
rm -rf nolib
rm -rf sharedlib
rm -rf staticlib

# use library mkmf.template
cp "$testdir"/mkmf.template.gfortran "$DART"/build_templates/mkmf.template

# regular filter build - this goes first so test directories get *.nc and obs_seq.out input files
cd "$DART"/models/lorenz_96/work
./quickbuild.sh clean
./quickbuild.sh mpif08 filter
rm *.o *.mod

#-----------------
# setup test direcories
cd "$testdir"
cp -r "$DART"/models/lorenz_96/work nolib
cp -r "$DART"/models/lorenz_96/work sharedlib
cp -r "$DART"/models/lorenz_96/work staticlib

#-----------------
# shared library build
cd "$DART"/developer_tests/library/shared/work
./quickbuild.sh clean
./quickbuild.sh mpif08

mpif90 -o filter "$DART"/assimilation_code/programs/filter/filter.f90 -I.  -L. -ldart

#-----------------
# static library build
cd "$DART"/developer_tests/library/static/work
./quickbuild.sh clean
./quickbuild.sh mpif08

mpif90 "$DART"/assimilation_code/programs/filter/filter.f90  -I. -L. -ldart  -L/$NETCDF_LIB -lnetcdff -o filter

#-----------------
# copy executables built from libraries to test directories
cd "$testdir"/sharedlib
rm filter
cp "$DART"/developer_tests/library/shared/work/filter .
cp "$DART"/developer_tests/library/shared/work/libdart.so .

cd "$testdir"/staticlib
rm filter
cp "$DART"/developer_tests/library/static/work/filter .

#-----------------
# run filter
cd "$testdir"/nolib/
mpirun -n 4 ./filter

cd "$testdir"/sharedlib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd) 
mpirun -n 4 ./filter

cd "$testdir"/staticlib/
mpirun -n 4 ./filter

check_bitwise "$testdir"/nolib/ "$testdir"/sharedlib/
check_bitwise "$testdir"/nolib/ "$testdir"/staticlib/

}

#-----------------
check_bitwise () {
diff -s "$1"/obs_seq.final "$2"/obs_seq.final

netcdffiles=(\
analysis.nc
filter_input.nc
filter_output.nc
perfect_input.nc
preassim.nc \
)

for f in "${netcdffiles[@]}"
do
  echo -n "$f" " "  
  nccmp -d "$1"/"$f" "$2"/"$f"
  echo "" 
done
}

main "$@"

