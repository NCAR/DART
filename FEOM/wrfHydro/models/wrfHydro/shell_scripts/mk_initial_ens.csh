#!/bin/csh

# this takes input HYDRO restart files from inDir for a given 
# yyyy mm dd hh
# and then generates nEns copies in smcPath (not user defined)
# Then it links the given fulldom file there and then runs the 
# fortran code to perturb the states in the files. You have to 
# edit the smcPerturb.f90 to get the right number of ensemble 
# members and other settings. 

## arguments: year, month, day, nEns, inDir, outDir, fulldom
## input files are 
# 1: 

if ($#argv != 7) then
    echo "Usage: yyyy mm dd nEns inDir outDir fulldom"
    exit 1
endif

set yy=$1
set mm=$2
set dd=$3
set nEns=$4
set restartPath=$5
set outPath=$6
set fulldom=$7

echo "Date: $yy $mm $dd"
echo "# of ensemble members: $nEns"
echo "Restart location: $restartPath"
echo "Output path: $outPath"
echo "Fulldom file: $fulldom"

cd $outPath || exit 1
pwd
rm -f restart*.nc

## right now we are only perturbing the hydro restart file.
set smcPath=/home/jamesmcc/fortranTools/SCRF/smcIcs/inout
rm -f ${smcPath}/* 

## make 10 copies, this has to correspond with whats in smcPerturb.f90
    foreach i (`seq 1 $nEns`)
    set ii=`printf "%04d" $i`
    \cp -v ${restartPath}/HYDRO_RST.${yy}-${mm}-{$dd}_00:00_DOMAIN* ${smcPath}/restart.hydro.${ii}.nc || exit 3
    \cp -v ${restartPath}/RESTART.${yy}${mm}${dd}00_DOMAIN* ${smcPath}/restart.${ii}.nc || exit 4
end

rm /home/jamesmcc/fortranTools/SCRF/smcIcs/Fulldom.symlink.nc
\ln -vs $fulldom /home/jamesmcc/fortranTools/SCRF/smcIcs/Fulldom.symlink.nc || exit 5

/home/jamesmcc/fortranTools/SCRF/smcIcs/smcPerturb.exe

cp ${smcPath}/* ${outPath}/. || exit 6


exit 0
