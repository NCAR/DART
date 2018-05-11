#!/bin/csh
#===============================================================================
# This script runs an ensemble forward in time. 
# The script can use perturbed forcing/parameters?

# Notes / Todo
# 1. Runs have to start and end on midnight currently.
# 2. Tried to make the entire script a torque job but I think using ex in a 
#    non-interactive terminal causes issues/failure. Trying new approach where
#    a torque script is cat'd to a file which is then run. 

#---------------------------------------------------------------------
## One application:
# Because streamflow cannot simply be perturbed, 
# this script takes a perturbed inital state and runs if forward
# with(?) or without perturbed forcing to generate an 
# ensemble of streamflow at an initial time.

#===============================================================================
# Calling directory:
#
# Input:
# Input arguments: endYyyy endMm endDD endHH inDir
# the inputDir should be of the form initialEnsemble.yyyymmdd.scheme 
# so that the start time can be grabbed from it.
#
# Output:
#===============================================================================

# Parse the input arguments
if ($#argv != 4 & $#argv != 5) then
    echo "Usage: endYyyy endMm endDD inDir <ppn>"
    exit 1
endif
set endYyyy = $1
set endMm   = $2
set endDd   = $3
set inDir   = $4
if ($#argv == 5) set ppn = $5
if(! $?ppn) set ppn = 16
if($ppn > 16) then 
    echo "ppn > 16, not allowed. Setting ppn=16."
    set ppn = 16
endif 

echo endYyyy: $endYyyy
echo endMm: $endMm
echo endDd: $endDd
echo inDir: $inDir

set wd = `pwd`
echo wd: $wd

# parse the inDir for the date
set inDirDate = `echo $inDir | cut -d. -f3`
set startYyyy = `echo $inDirDate | cut -c1-4`
set startMm = `echo $inDirDate | cut -c5-6`
set startDd = `echo $inDirDate | cut -c7-8`

echo $startYyyy

set startDateSec = `date -d $startYyyy-$startMm-$startDd +%s`
set endDateSec = `date -d $endYyyy-$endMm-$endDd +%s`
@ newKhour = ($endDateSec - $startDateSec) / 3600

set outPath = initialEnsembleWFlow.$endYyyy$endMm$endDd.`echo $inDir | cut -d. -f3`

set num_states = `ls $inDir/restart.hydro.*.nc -1 | wc -l`


echo "Start Date: $startYyyy $startMm $startDd"
echo "End Date: $endYyyy $endMm $endDd"
echo "Number of hours to advance: $newKhour"
echo "Output path: $outPath"
echo "Number of ensemble members: $num_states"
# the following follows advance_model.csh to some degree.

# mv restarts from output to here
\mkdir -p $outPath || exit 1  
# but actually run here
\mkdir -p RUNS.$outPath || exit 1
cd RUNS.$outPath

## get the name lists
\cp ../namelist.hrldas .
\cp ../hydro.namelist .

# Get parameters.
foreach FILE ( ../PARAMS.gathered/* ) 
    echo $FILE
#   \cp -v ../$FILE . || exit 2 ## these dont change so link them instead.
    \ln -sf $FILE . || exit 2  
end

## write the torque script to run the model in parallel
echo "#\!/bin/csh" >! jobScript.csh
echo "#PBS -l nodes=1:ppn=$ppn,walltime=100:01:00" >> jobScript.csh
echo "#PBS -e stderr.txt" >> jobScript.csh
echo "#PBS -o stdout.txt" >> jobScript.csh
echo "#PBS -N ini_ens_fwd" >> jobScript.csh
echo "source ~/.cshrc" >> jobScript.csh
echo "cd `pwd`" >> jobScript.csh
echo "\rm -f qsubIsDone.criticalFile" >> jobScript.csh
echo "mpirun -np $ppn ../wrf_hydro.exe >& /tmp/iniEnsFwdJunk" >> jobScript.csh
echo "\rm /tmp/iniEnsFwdJunk" >> jobScript.csh
echo "touch  qsubIsDone.criticalFile" >> jobScript.csh
echo "exit 0" >> jobScript.csh

## generate the random mean of each perturbation
echo "../randomSeq.rsh uniform $num_states .5 1.5"
\cp ../randomSeq.rsh .
set randomMeans = `./randomSeq.rsh uniform $num_states .5 2`
echo The random mean precip multipliers:
echo $randomMeans

# Loop through each state /ensemble member
set state_copy = 1
while($state_copy <= $num_states)

    set instance        = `printf "%04d" $state_copy`
    echo "ensemble member: $instance"

    \rm -f  restart.nc  restart.hydro.nc 
    \rm -f  HYDRO_RST.*  RESTART.* 
    \rm -f  *.LDASOUT_DOMAIN*
    \rm -f  *.LSMOUT_DOMAIN*  *.RTOUT_DOMAIN*  *.CHRTOUT*  *.CHANOBS*  frxst_pts_out.txt
    \rm -f  qstrmvol*  diag_hydro.*  stderr.txt stdout.txt  GW_*.txt  *.GW_DOMAIN* 

    # need the wrfHydro restart files 
    \ln -sv ../$inDir/restart.$instance.nc  restart.nc   || exit 2
    \ln -sv ../$inDir/restart.hydro.$instance.nc  restart.hydro.nc   || exit 2

    # time stuff
    # KHOUR
    set numadvances = $newKhour

    # Have to match the start time in hrldas with the restart files,
    # else the forcing data seems to have no effect.
    set restartFileTime = `ncdump -v Times restart.nc | tail -2 | head -1 | cut -d'"' -f2`
    set restartFileYyyy = `echo $restartFileTime | cut -d- -f1`
    set restartFileMm = `echo $restartFileTime | cut -d- -f2`
    set restartFileDd = `echo $restartFileTime | cut -d- -f3 | cut -d_ -f1`
    set restartFileHh = `echo $restartFileTime | cut -d_ -f2 | cut -d: -f1`

## also need to check for /change  output dir for noah??

ex namelist.hrldas <<ex_end
g;KHOUR ;s;= .*;= $numadvances;
g;START_YEAR;s;= .*;= $restartFileYyyy;
g;START_MONTH;s;= .*;= $restartFileMm;
g;START_DAY;s;= .*;= $restartFileDd;
g;START_HOUR;s;= .*;= $restartFileHh;
wq
ex_end


    ## Perturb forcing 
    ## (0. outside this loop (state_copy) define a uniform random sequence of mean precip adjustment)
    ## 1. Edit the path to the forcing dir in the namelist.hrldas.
ex namelist.hrldas <<ex_end
g;INDIR ;s;= .*;= "FORCING.perturbed.spinup/";
wq
ex_end

    ## 2. for all forcing files in the time period
    @ firstForcingFileInd = `\ls -1 ../FORCING/ | \grep -n $startYyyy$startMm${startDd}00 | \
                             cut -d: -f1`
    set forcingFiles = `\ls -1 ../FORCING/ | \tail -n+$firstForcingFileInd | \head -$numadvances`
    ## 3. use the mean precip adjustment to define a sequence of normal errors about that mean
    ## through time. (use randomSeq.rsh, which computes text so you dont have to do it in shell)
    set timePerts = `../randomSeq.rsh uniform $numadvances \
			"$randomMeans[$state_copy]*(1-.1)" \
                        "$randomMeans[$state_copy]*(1+.1)"`

    @ countForc=1
    \mkdir FORCING.perturbed.spinup
    foreach fForcing ( $forcingFiles )
	##  multiply the precip in the copy using ncap2 and write a  
	##     copy to ../FORCING.perturbed.spinup.iEns.endYyyyendMmendDd 
	set iPerturb = $timePerts[$countForc]

	ncap2 -s "RAINRATE=RAINRATE*${iPerturb}" ../FORCING/$fForcing \
		 FORCING.perturbed.spinup/$fForcing

	@ countForc++
    end 

\ls

    echo "advance the model"
	qsub jobScript.csh 
    echo "waiting for qsub model run to finish"
    while ( ! -e qsubIsDone.criticalFile )
	sleep 1
    end
    echo "qsub model run completed"    
    #mpirun -np $ppn ./wrf_hydro.exe >& /tmp/iniEnsFwdJunk
    #../wrf_hydro.exe >& /tmp/jamesmccWfrHydroEnsOutputJunk
    #\rm -f /tmp/jamesmccWfrHydroEnsOutputJunk
    
    ## keep the desired restart files
    \mv -v `ls RESTART* | tail -1` ../$outPath/restart.$instance.nc || exit 4
    \mv -v `ls HYDRO_RST* | tail -1` ../$outPath/restart.hydro.$instance.nc || exit 5

    ## remove the undesired ouput files
    \rm -rf *.LDASOUT_DOMAIN* *.CHANOBS_DOMAIN* *.CHRTOUT_DOMAIN* \
	   *.LSMOUT_DOMAIN1 *.RTOUT_DOMAIN1 diag_hydro.* \
	   frxst_pts_out.txt GW_*.txt qstrmvolrt_accum.txt \
	   RESTART.*_DOMAIN* HYDRO_RST.*_DOMAIN* qsubIsDone.criticalFile \
	   stderr.txt stdout.txt

    ## clean up the forcing
    \rm -rf FORCING.perturbed.spinup

\ls

    @ state_copy++
end

cd ..


exit 0

