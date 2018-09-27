#!/bin/bash

## With no argument, this script runs ALL tests for the WRF Hydro model.
## With integer argument, this script runs that specified test. 

## Setup which test cases will be performed.
nArgs="$#"

testCases=(0)

if (( $nArgs > 0 )); then
    testCases=("$@")
fi

if (( ${testCases[0]} == 0 )); then
    testCases=(`ls -d test* | egrep 'test[1-9](1-9)?$'`)
else 
    for tt in `seq 0 $((nArgs-1))`

    do
        testCases[$tt]='test'${testCases[$tt]}
    done
fi

##=============================================================================================
function doTest {
    ## This function performs the test

    ## check that we are in models/wrfHydro/work before doing anything.
    curDir=`pwd`
    inWork=`echo $curDir | grep 'models/wrfHydro/work'`
    if [[ ! "$inWork" ]]; then 
        return 1
    fi

    ## parse the single arg
    nReqArgs=1
    if [ "$#" -ne $nReqArgs ]; then
        echo 'doTest requires '$nReqArgs' arguments, you supplied '"$#"'.'
        return 1
    fi
    testDir=$1
    
    echo "||||||||||||||===================================================||||||||||||||"
    echo `date`

    ## TESTS ARE RUN IN THE INDIVIDUAL TEST DIRS
    ## Each test dir is required to have 
    ## input.nml
    ## obs_seq.out
    ## initialEnsemble directory
    ## DOMAIN directory
    ## FORCING directory
    ## PARAMS.noah or PARAMS.noahMP directories
    ## validation/Prior_Diag.nc validation/Posterior_Diag.nc

    ## The following need placed up one directory
    runReqs=( DOMAIN FORCING PARAMS.noah PARAMS.noahMP PARAMS.hydro )
    for rr in "${runReqs[@]}"
    do
        echo $rr
        rm -f $rr
        ln -sf $testDir/$rr . || return 1
    done

    cd $testDir || return 1
    exes=( setup_filter.csh wrf_hydro.exe )
    for ee in "${exes[@]}"
    do
        echo $ee
        ln -sf ../$ee . || return 1
    done

    ## testDescription.txt
    ## cat description to log

    echo "./setup_filter.csh"
    ./setup_filter.csh || return 1
    
    echo "./run_filter.csh"
    ./run_filter.csh || return 1

    ## how to verfiy in a standard way?
    ## each test directory has it's own validate.sh which only echos output if any sub test fails. 
    echo "./validate.sh"
    validOut=`./validate.sh`

    if [[ ! $validOut == '' ]]
    then
        echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' | tee /dev/stderr 
        echo "Test $testDir fails. Output from $testDir/validate.sh:"  | tee /dev/stderr
        echo $validOut                                                 | tee /dev/stderr 
        echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' | tee /dev/stderr 
        echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' | tee /dev/stderr 
        return 1
    fi

    echo "$testDir PASSED the tests" | tee /dev/stderr
    echo 

    cd ../
    
    return 0
}

##=============================================================================================

for tt in "${testCases[@]}"
do 
    echo $tt logged to ${tt}.log
    ## only send stderr to the terminal, log both stdout and stderr
    doTest $tt > ${tt}.log 2> >(tee ${tt}.log >&2)
done

exit 0
