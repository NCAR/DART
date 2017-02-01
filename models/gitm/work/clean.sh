#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# From Alexey Morozov

# This script assumes we are in the work dirctory.
# It spins up noe=$1 (1st argument to the script) ensemble members +
# (if ts=$3=1) one truth simulation with f1t=$4 f107 value used for
# truth (if ts=1) and f107 sampled from a normal distribution with
# mean f1m=$5 and standard deviation f1s=$6 for the ensemble members

# ./clean.sh 1 1 0 150.0 170.0 1.0 

export noe=$1 # number of DART ensemble members
export nop=$2 # number of processors per each member
export  ts=$3 # need to prespin truth? 1=yes 0=no
export f1t=$4 # mean of the f107 initial ensemble distribution
export f1m=$5 # mean of the f107 initial ensemble distribution
export f1s=$6 # standard deviation (sigma) of the f107 initial ensemble distribution

rm *~                # all the emacs backups
rm -rf advance_temp* # the old folders

rm perfect_ics       # perfect_ics will be created by this script if ts=1
rm filter_ics*       # filter_ics  will be created by this script
rm filter_control*   # old DART control files
rm filter_lock*      # old DART lock files

# create a file (fens.txt) containing the initial distribution 
# of f107 for the ensemble members
python ../python/create_init_distr.py $noe $f1m $f1s fens.txt

# create the temporary directories -
# one per ensemble member and start each prespin run in the background
for (( i = 1 ; i <= $[$noe + $ts]; i++ ))
do
    curdir='advance_temp_e'$i #name
    echo "clean.sh: right now we are creating folder" $curdir

    #the following is a modified version of GITM's Makefile's "make rundir"
    mkdir -p ${curdir}/UA
    cd ${curdir}
       ln -s ../../GITM2/util/EMPIRICAL/srcIE/data EIE
       ln -s ../../GITM2/src/PostProcess.exe PostGITM.exe
       cd UA
	   mkdir restartOUT restartOUT.out data
	   ln -s restartOUT.out restartIN
	   ln -s ../../../GITM2/src/pGITM .
 	   ln -s ../../../GITM2/srcData DataIn
           # to save on inodes (number of files) quota on Pleiades,
           # instead of linking every single file in this dir we link just the folder
           # mkdir DataIn #not needed (compare to the "make rundir" section in the Makefile)
           # cd DataIn
           #    ln -s ../../../../GITM2/srcData/* .
           #    rm -f CVS
           # cd ..
       cd ..
       chmod -R 755 UA || exit 1
       ln -s ../../GITM2/src/GITM.exe .
       ln -s UA/* .

       if [ $i -eq $[$noe + 1] ] ; then #truth folder
	
           # the data part of UAM.in - how often to write output files etc.
           # 1) (this truncated version will be augmented with more lines below to become a "full" UAM.in)
           # 2) the time part of UAM.in telling GITM for how long to prespin
           # 3) gitm_to_dart and dart_to_gitm look for UAM.in there as well but only to know the number of blocks in lat and lon
           # 4) this script assumes that all the dat files (satellites, imf, power, etc) are in the
           #    DART/models/gitm/GITM2/run folder and copies them from there to all the ensemble folders.
	   cp ../../GITM2/srcData/UAM.in.test.DART.data_truth   UAM.in.truncated
	   cp ../../GITM2/srcData/UAM.in.test.DART.time_prespin .
	   cp UAM.in.truncated UA/restartOUT/UAM.in
	   cp ../../GITM2/run/*.dat .

           echo "#DART"   >> UAM.in.truncated
           # this is just a plain GITM run (to collect truth data or just see how bad the bias is)
           echo "0  useDART (0=no, 1=master, 2=slave)"  >> UAM.in.truncated
           echo " "       >> UAM.in.truncated
           echo "#F107"   >> UAM.in.truncated
           echo $f1t      >> UAM.in.truncated #what fixed f107 you want for this truth run?
           echo $f1t      >> UAM.in.truncated #what fixed f107a (81day centered average) you want?
           # If you don't want f107 and f107a fixed, comment the above 3 lines and uncomment the next 2 lines
	   # echo "#NGDC_INDICES"      >> UAM.in.truncated
           # echo "UA/DataIn/f107.txt" >> UAM.in.truncated #use time varying f107 from this text file instead
           echo " "       >> UAM.in.truncated

           # just write to the screen that the truth folder just got created
	   echo "clean.sh: truth folder got created"

       else #ensemble members (not truth simulation) folders

           # 1) copy the "truncated" version of UAM.in
           #    (this truncated version will be augmented with more lines below to become a "full" UAM.in)
           # 2) the time part of UAM.in telling GITM for how long to prespin
           # 3) DART wants a UAM.in there too to read the grid data (NBlocksLon/Lat)
           # 4) this script assumes that all the dat files (satellites, imf, power, etc)
           #    are in the DART/models/gitm/GITM2/run folder and copies them from there to all the ensemble folders.
	   cp ../../GITM2/srcData/UAM.in.test.DART.data_ens UAM.in.truncated
	   cp ../../GITM2/srcData/UAM.in.test.DART.time_prespin .
	   cp UAM.in.truncated UA/restartOUT/UAM.in
	   cp ../../GITM2/run/*.dat .

           echo "#DART"   >> UAM.in.truncated
           # for prespin purpose, everybody is a "true simulation" in a sense that nobody needs to
           # be inflated - just keep your initial f107 value untill we start assimilating
           echo "0  useDART (0=no, 1=master, 2=slave)"  >> UAM.in.truncated
           echo " "       >> UAM.in.truncated
           # fens.txt got created by the python call to ../../create_init_distr.py above,
           # here just read the i-th line of that file
	   export f1cur=$(sed -n $i'p' ../fens.txt)
           echo "#F107"   >> UAM.in.truncated
           echo $f1cur    >> UAM.in.truncated #f107 for this day
           echo $f1cur    >> UAM.in.truncated #f107a - average over 81 days
           echo " "       >> UAM.in.truncated

           # just write to the screen what ensemble just got processed
           echo "clean.sh: ensemble" $i "has f107 of" $f1cur
	
       fi

       cat UAM.in.test.DART.time_prespin >  UAM.in #put the time part of UAM.in into UAM.in
       cat UAM.in.truncated >> UAM.in #put the extended data part (with #DART and #F107) of UAM.in into UAM.in
       echo "#RESTART"      >> UAM.in
       echo "F"             >> UAM.in #prespin is the first time we ran this GITM ever, so it can't be a restart.
       echo " "             >> UAM.in
       echo "#END"          >> UAM.in

       #Add restart is true to the truncated version for later use
       echo "#RESTART"      >> UAM.in.truncated
       echo "T"             >> UAM.in.truncated #after prespin is done, this can be restarted
       echo " "             >> UAM.in.truncated

       # copy DART programs into each dir
       ln -s ../dart_to_gitm . || exit 2 # converts DART files into GITM files
       ln -s ../gitm_to_dart . || exit 3 # the inverse
       cp ../input.nml .       || exit 4 # control info for dart_to_gitm and gitm_to_dart

       # we need to change gitm_restart_dirname, as we moved down one directory.
       # remove "advance_temp_e1/" from gitm_restart_dirname
       sed -i'' 's#advance_temp_e1/##' input.nml

       if [ $?LSB_JOBID ] ; then
          echo "Runnin under LSF in "`pwd`
          mpirun.lsf ./GITM.exe &
       else
          # "hf" is used to tell the mpiexec to run on specific (separate) nodes
          # start each ensemble member in the background and save the stdout and stderr
          # into a file called ens_ps.txt in each folder
          head -n $[$nop*$i] $PBS_NODEFILE | tail -n $nop > hf
          echo "Runnin under PBS on nodes ... "
          cat hf
          (mpiexec -hostfile hf -n $nop GITM.exe >& ens_ps.txt ) &
       fi

    cd ..

done

wait

#now iterate over all ensemble members with the goal of seeing if they all finished
# prespinning already and then convert GITM files to DART files and do some cleaning up
for (( i = 1 ; i <= $[$noe + $ts]; i++ ))
do

 curdir='advance_temp_e'$i
 cd  $curdir  || exit 5

    test -e GITM.DONE
    export s=$?
    while [ $s -eq 1 ] #while it doesn't exist
    do
        echo "clean.sh: GITM.DONE doesn't exist in " $curdir
        sleep 5
        test -e GITM.DONE
        export s=$?
    done

    echo "listing contents of directory $curdir"
    ls -lR

    rm GITM.DONE               # so that adv_model doesn't trip on old G.D files
    rm UAM.in                  # remove the used full UAM.in
    cp UAM.in.truncated UAM.in # prepare a truncated UAM.in for the next assimilation window

    rm UA/data/3D*.*           # so that there is less stuff to move back and forth
    rm UA/data/iriOut_*.*      # so that there is less stuff to move back and forth
    rm UA/data/log*.*          # so that there is less stuff to move back and forth

    ./gitm_to_dart || exit 6   # convert GITM restart files to DART initial file

    echo "listing contents of directory after gitm_to_dart $curdir"
    ls -lrt

    if [ $i -eq $[$noe + 1] ] ; then #perfect IC
	mv -v dart_ics ../perfect_ics
        # this file is used by ./perfect_model_obs, but I don't use it to get true data -
        # I just run GITM another time - see pbs_file.sh
    else  #filter ICs
	if [ $i -lt 10 ] ; then
            # filter ICs with index less than 10
            # these files are used by ./filter
	    mv -v dart_ics ../filter_ics.000$i
	else  #f_ICs with index greater or equal to 10 but less than 99
	    mv -v dart_ics ../filter_ics.00$i
	fi
    fi

 cd ..

done

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

