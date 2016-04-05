#!/bin/csh

# Because streamflow cannot simply be perturbed, 
# this script takes a perturbed inital state and runs if forward
# with(?) or without perturbed forcing to generate an 
# ensemble of streamflow at an initial time.

# input arguments endYyyy endMm endDD endHH inDir
# the inputDir should be of the form initialEnsemble.yyyymmdd.scheme 
# so that the start time can be grabbed from it.

## note, configure for both qsub and without (see run_filter.csh)

#-----------------------------------------------------------
if ($#argv != 6 & $#argv != 7) then
    echo "Usage: endYyyy endMm endDD endHH inDir(relative to calling dir) ppn inputNml"
    exit 1
endif
set endYyyy = $1
set endMm   = $2
set endDd   = $3
set endHh   = $4
set inDir   = $5
set ppn     = $6
if ($ppn > 16) then 
    echo "Right now only configured to use a single node, setting ppn=8."
    set ppn = 8
endif 
if ($#argv == 7) set inputNml = `readlink -f $7`


#-----------------------------------------------------------
## Examine the contents of the inDir and its internal consistency
## The number of restart files.
set nLsmRestarts       = \
    `ls -1 $inDir/restart* | egrep 'restart.[0-9]*.nc' | wc -l`
set nHydroRestarts     = \
    `ls -1 $inDir/restart.hydro* | egrep 'restart.hydro.[0-9]*.nc' | wc -l`
set nAssimOnlyRestarts = \
    `ls -1 $inDir/restart.assimOnly* | egrep 'restart.assimOnly.[0-9]*.nc' | wc -l`
@ assimOnlyActive = ( $nAssimOnlyRestarts != 0 )
if ( $nLsmRestarts != $nHydroRestarts ) then 
    echo "The number of LSM restart files does not match the number of hydro restart files."
    exit 1
endif

if ( $assimOnlyActive ) then 
    if ( $nLsmRestarts != $nAssimOnlyRestarts ) then
	echo "The number of LSM restart files does not match the number of assimOnly restart files."
	exit 1
    endif
    if ( ! $?inputNml ) then
        echo "restart.assimOnly.nc files detected, but no input.nml file has been specified."
        exit 1
    endif
endif

set nEns = $nLsmRestarts

## Check that they *all* have the same restart times.
set restartTime = `ncdump -v Times $inDir/restart.0001.nc | tail -2 | head -1 | cut -d'"' -f2`
foreach iEns ( `seq 1 $nEns` )
    set iEnsFmt      = `printf "%04d" $iEns`
    set rstTimeLsm = \
	`ncdump -v Times $inDir/restart.${iEnsFmt}.nc | tail -2 | head -1 | cut -d'"' -f2`
    set rstTimeHydro = \
	`ncdump -h $inDir/restart.hydro.${iEnsFmt}.nc | grep Restart_ | cut -d'=' -f2 | tr -d ' ";'`
    if ( $assimOnlyActive ) then 
    set rstTimeAssimOnly = \
	`ncdump -h $inDir/restart.assimOnly.${iEnsFmt}.nc | grep ':Restart_' | cut -d'=' -f2 | tr -d ' ";'`
    endif 
    if ( $restartTime != $rstTimeLsm       |  \
         $restartTime != $rstTimeHydro     ) then 
	 echo "There are non-matching restart times between the lsm or hydro restart files."
	 exit 2
    endif


    if ( $assimOnlyActive ) then 
	if ( $restartTime != $rstTimeAssimOnly ) then 
	    echo "There are non-matching restart times between the lsm and assimOnly restart files."
	    exit 2
	endif 
    endif
end

# Set the restart time
set startYyyy = `echo $restartTime | cut -d- -f1`
set startMm   = `echo $restartTime | cut -d- -f2`
set startDd   = `echo $restartTime | cut -d- -f3 | cut -d_ -f1`
set startHh   = `echo $restartTime | cut -d_ -f2 | cut -d: -f1`

# Calculate the khour from the restart time and supplied end time
set startDateSec = `date -ud "UTC $startYyyy-$startMm-$startDd $startHh hour" +%s`
set endDateSec = `date -ud "UTC $endYyyy-$endMm-$endDd $endHh hour" +%s`
@ newKhour = ($endDateSec - $startDateSec) / 3600

# The output path
set outPath = ensembleRun.$startYyyy$startMm$startDd-$endYyyy$endMm$endDd

# All setup. Print the details. 
echo "Start Date: $startYyyy $startMm $startDd $startHh"
echo "End Date  : $endYyyy $endMm $endDd $endHh"
echo "Number of hours to advance: $newKhour"
echo "Number of ensemble members: $nEns"
echo "assimOnlyActive: $assimOnlyActive"
echo "Input dir     : $inDir"
echo "Output/Run dir: $outPath"
echo "ppn : $ppn"

# Run ensembles within this dir
if ( -e $outPath ) then 
    echo 'Output/Run directory already exists, exiting'
    exit 1
endif
\mkdir -p $outPath || exit 1
cd $outPath
ln -sf ../$inDir .  ## make this symlink so it's clear what these runs started from. 

## get the name lists
\cp ../namelist.hrldas .
\cp ../hydro.namelist .
\ln -sf ../DOMAIN .

# Have to match the start time in hrldas with the restart files,
# else the forcing data seems to have no effect.
ex namelist.hrldas <<ex_end
g;KHOUR ;s;= .*;= $newKhour;
g;START_YEAR;s;= .*;= $startYyyy;
g;START_MONTH;s;= .*;= $startMm;
g;START_DAY;s;= .*;= $startDd;
g;START_HOUR;s;= .*;= $startHh;
wq
ex_end
# Note that the namelists used for dat assume the restart files are
# named restart.nc and restart.hydro.nc, so that specification in the
# namelist dosent need need edited to reflect the restart time.

## get the forcing files before changing dirs
set forcDir = `grep INDIR namelist.hrldas | tr -d ' "' | tr -d "'" | grep -v '^[\!\]' | cut -d'=' -f2`
if ( `echo $forcDir | egrep '^[.]' | wc -l` ) set forcDir = `echo ../${forcDir}`
ln -s $forcDir .

\cp ../setup_ensemble_run.csh .

#===============================================================================
# Loop through the ensemble members
set iEns = 1
while($iEns <= $nEns)

    set iEnsFmt        = `printf "%04d" $iEns`
    echo "ensemble member: $iEnsFmt"

    ## the rest of the insde of this loop could be put in a separate file to || it.
    ## can save time when many forcing files are to be adjusted
    if ( `ls -1 ensemble.*/wrfHydroStillWorking.dum | wc -l` >= $ppn ) then
	echo ' ********** Waiting for processors... ********** '
	echo ' ********** Waiting for processors... ********** '
	echo ' ********** Waiting for processors... ********** '
    endif
    while ( `ls -1 ensemble.*/wrfHydroStillWorking.dum | wc -l` >= $ppn )
	sleep 1
    end
   
    set ensDir = ensemble.${iEnsFmt}
    \mkdir $ensDir
    ## this control file gets cleaned up in run_wrfHydro.csh
    touch $ensDir/wrfHydroStillWorking.dum

    ./setup_ensemble_run.csh $inDir $iEnsFmt $assimOnlyActive $inputNml &

    @ iEns++
end

## we dont want this script to finish until all the runs are complete. 
while ( `ls -1 ensemble.*/wrfHydroStillWorking.dum | wc -l` > 0  )
    sleep 1
end



exit 0

