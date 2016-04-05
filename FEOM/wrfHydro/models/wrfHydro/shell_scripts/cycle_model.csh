#! /bin/csh
## cycling is done relative to a previous assimilation run.
## It relies on an OUTPUT directory from that  with restart files
## If "precycling" is desired, additonal restart files are required. Precycling 
## is running the model from restarts which existed prior to the assimilation for 
## times at the beginning of the assimilation run (assimilation may not be desired. I 
## could also configure this to do assimilation from earlier restarts, but to only use 
## precipitation multiplier at a give moment in the cycling.)

## Flow of Cycling
## (restarts are timestamped with the time to be run)
## Required Inputs: 
##   1: nStepsCycle - if nStepsCycle=3, then after running t=3(->4), go back and re run t=1,2,3.
##      (given nStepsCycle, the following restart timestamps are true:
##         after running time t+nStepsCycle-1(->t+nStepsCycle), 
##           run timesteps, t(->t+1), ... ,t+nStepsCycle-1(->t+nStepsCyclen) 
##       )
##   2: prevOutputDir/ - which contains the following:
##        + initial_restarts/ :restart.iEns.nc@t=1, restart.hydro.iEns.nc@t=1
##        + model_integration.{t}-{t+1}.iEns/
##           +restart.assimOnly.iEns.nc@t+1
##   3: precycleModelRestarts/ ONLY IF precycling, then also need to specifiy this dir 
##      containing model (not assimOnly) restarts prior to the assimilation output in
##       prevOutputDir. 

## Flow:
## say initial_restarts/ :restart.iEns.nc@t=1, restart.hydro.iEns.nc@t=1 (= run from t=0(->1))
## The first restart.assimOnly.iEns.nc is for t=1, prior to model advance in assimilation.
## This file should be in this directory.
## (Here 1 is the first model integration performed during assimilation and the timestamp of 
## the input restart files).
## PRECYCLING: precyclemodelRestarts/       restart.iEns.nc@t=1-nStepsCycle, 
##                                    restart.hydro.iEns.nc@t=1-nStepsCycle

## ***
## fix assim put the restart.assimOnly  in t+1 instead of t. delete last and fix.
## also need to fix (??) in advance model.
## for consistency, the restart.assimOnly.nc file should be associated with the timestep before 
## the time it "forces". The first one should go in initial_restarts.
## ***


echo $nStepsCycleIn
echo $prevOutputDir

## these are the times given to namelist.hrldas
set nTimes = `ls -1d $prevOutputDir/model_integration.*001 |wc -l`
set cycleEnd = `seq 1 $nTimes`
set nRewind = `echo $cycleEnd`
set cycleStart = `echo $cycleEnd`

## set cycleEnd times
set count = 0
foreach tt ( `ls -1d $prevOutputDir/model_integration.*001` )
    @ count++
    set cycleEnd[$count] = `echo $tt | cut -d. -f5 | cut -d- -f1`
end

## set cycleStart times
set count = 0
foreach tt ( `seq 1 $nTimes` )
    @ count++
    set nRewind[$count] = `echo $count $nStepsCycleIn | xargs -n 1 | sort -n | head -1`
    @ startInd = $count - $nRewind[$count] + 1
    set cycleStart[$count] = $cycleEnd[$startInd]
end

## reset the symlink for the initial ensemble
\rm -rf initialEnsemble
\rm -rf initialEnsemble.CYCLE
mkdir initialEnsemble.CYCLE
\ln -sf initialEnsemble.CYCLE initialEnsemble
\cp precycleAssimOnlyInit/* initialEnsemble/

\rm Prior_Diag.CYCLE.nc Posterior_Diag.CYCLE.nc

## cycle
foreach tt ( `seq 1 $nTimes` )
    
    echo '==================================================================='
    echo

    set startY = `echo $cycleStart[$tt] | cut -c1-4`
    set startM = `echo $cycleStart[$tt] | cut -c5-6`
    set startD = `echo $cycleStart[$tt] | cut -c7-8`
    set startH = `echo $cycleStart[$tt] | cut -c9-10`

    set endY = `echo $cycleEnd[$tt] | cut -c1-4`
    set endM = `echo $cycleEnd[$tt] | cut -c5-6`
    set endD = `echo $cycleEnd[$tt] | cut -c7-8`
    set endH = `echo $cycleEnd[$tt] | cut -c9-10`

    echo First : Last
    echo $startY $startM $startD $startH
    echo $endY $endM $endD $endH

    ## move the assimOnly files previous obtained for this timestep into initialEnsemble
    ## this is confusing b/c it's CYCLING you actually want to use the restart at time+1
    set currentAssimOnlies = \
	    `ls -d $prevOutputDir/model_integration.*$endY$endM$endD$endH.*/restart.assimOnly.nc`
    echo currentAssimOnlies $currentAssimOnlies
    foreach ii ( $currentAssimOnlies )
	set iens=`echo $ii | cut -d. -f6 | cut -d/ -f 1`
	cp $ii initialEnsemble/restart.assimOnly.$iens.nc
    end

    ## run the ensemble forward without assim to get new restarts??
    ## right now this is hardcoded
    if ( $nRewind[$tt] < $nStepsCycleIn ) then
	echo --- PRECYCLE ---
	echo --- First : Last ---
	@ startHPre = $endH - $nStepsCycleIn + 1
	echo $startY $startM $startD $startHPre
	echo $endY $endM $endD $endH
	./precycle_ens.csh $startY $startM $startD $startH \
	    initialEnsemble.$startY$startM$startD$startHPre.toPreCycle/ precycleAssimOnly/ 1

	# move the restarts to the initial restarts
	# again hardcoded
	set preRunsDir = RUNS.initialEnsembleWFlow.$endY$endM$endD$endH.toPreCycle/
	echo $preRunsDir
	foreach ff ( $preRunsDir/ensemble.*/RESTART*$endY$endM$endD${endH}* )
	    set iens =  `echo $ff | cut -d. -f5 | cut -d/ -f1` 
	    cp $ff initialEnsemble/restart.$iens.nc
	end 

	foreach ff ( $preRunsDir/ensemble.*/HYDRO_RST*$endY-$endM-${endD}_$endH* )
	    set iens =  `echo $ff | cut -d. -f5 | cut -d/ -f1` 
	    cp $ff initialEnsemble/restart.hydro.$iens.nc
	end 

	## the second? state generated should be put in to an "initialENsemble.date.toPrecycle" dir.

    else 
	## if not precycling just mv the working restarts to the initialEnsembles folder
        mv restart.0*.nc restart.h*.nc restart.a*.nc initialEnsemble/
    endif 

## need to setup obs_seq.out!!!!

    # assimilate
    ./setup_filter 
    ./run_filter

    # manage OUTPUT? I think the last output files written will be what we want. 
    ## hope for no conflict
    if (! -e Prior_Diag.CYCLE.nc ) 
	mv Prior_Diag.nc Prior_Diag.CYCLE.nc
	mv Posterior_Diag.nc Posterior_Diag.CYCLE.nc
    else 
	ncks -F -d -A time,$nRewind Prior_Diag.nc Prior_Diag.CYCLE.nc 
	ncks -F -d -A time,$nRewind Posterior_Diag.nc Posterior_Diag.CYCLE.nc 
    endif 

end 







