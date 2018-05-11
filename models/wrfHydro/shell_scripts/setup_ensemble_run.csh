#!/bin/csh

set inDir           = $1
set iEnsFmt         = $2
set assimOnlyActive = $3
set inputNml        = $4

set ensDir = ensemble.${iEnsFmt}
cd $ensDir
## get the name lists
\cp ../namelist.hrldas .
\cp ../hydro.namelist .

# Get parameters.
foreach FILE ( ../../PARAMS.gathered/* ) 
	echo $FILE
	\ln -sf $FILE . || exit 2  
end


# the wrfHydro restart files for the individual ensemble members
\ln -sv ../$inDir/restart.$iEnsFmt.nc  restart.nc   || exit 2
\ln -sv ../$inDir/restart.hydro.$iEnsFmt.nc  restart.hydro.nc   || exit 2
if ( $assimOnlyActive ) then 
	ln -sv ../$inDir/restart.assimOnly.$iEnsFmt.nc  restart.assimOnly.nc   || exit 2
endif

## If assimOnly variables need applied, this script does it.
if ( $assimOnlyActive ) then 
        ln -sf $inputNml input.nml
	cp ../../apply_assimOnly.csh .
	./apply_assimOnly.csh
	set assimOnlySuccess = $?
	if ( $assimOnlySuccess != 0 ) then
	    echo "apply_assimOnly.csh failed with return value $assimOnlySuccess"
	    exit 3
	endif
endif 
   
## log 
cp ../../log_meta_wrfHydro.csh .
./log_meta_wrfHydro.csh

## Actually run
cp ../../run_wrfHydro.csh .
./run_wrfHydro.csh $iEnsFmt &


exit 0
