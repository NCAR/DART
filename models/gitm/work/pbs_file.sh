#!/bin/bash
#
# DART $Id$
#
#PBS -S /bin/bash
#PBS -N dart_test_tec
#PBS -l procs=21,pmem=1000mb,walltime=1:00:00
#PBS -A ACCOUNT
#PBS -l qos=ACCOUNT
#PBS -q QUE
#PBS -M USERNAME@EXAMPLE.com
#PBS -m abe
#PBS -V

#No semicolons after numerical values please, and a space before comments (this file is read into matlab, so certain rules have to be followed) also no / in the end of paths as I add them later
hd=~/HEAD/DART/models/gitm #directory where gitm is located
rd=/nobackup/NYXUSERNAME #temporary directory where you want to run the simulation on
noe=20 #number of ensemble members
nop=1 #number of processors per member

prs=1        #do you want to prespin? 1=yes, 0=no. Note, if you say no, it will assume you have "prespun" directories already available in folder specified on line 160ish below, so default, do Yes=1, even for short period of time just to get the folder structure in place (to adjust the length of prespin, change $hd copy of gitm/GITM2/srcData/UAM.in.test.DART.time_prespin ) 
ts=1         #do you want to do a "Truth Spin"? Yes = 1, no = 0. Not necessary for anything but plotting and generating "truth data" for later runs (so setting it to 0 is ok for initial tests).
est_f107=1   #do you want to estimate f107? 0=no, 1=yes. If you don't estimate f107, it is just kept constant at whatever it was initialized in clean.sh (will overwrite a line in the $hd of model_mod.f90 before building)
nbn=2        #Number of Blocks loNgitude (if bui=1, will overwrite a line in the $hd copy of gitm/GITM2/srcData/UAM.in.test.DART.data_ens and gitm/GITM2/srcData/UAM.in.test.DART.data_truth)
nbt=2        #Number of Blocks laTitude (see line before) 
ncn=9        #Number of Cells in each block in loNgitude (if bui=1, will overwrite a line in the $hd copy of gitm/GITM2/src/ModSize.f90)
nct=9        #Number of Cells in each block in laTitude (if bui=1, will overwrite a line in the $hd copy of gitm/GITM2/src/ModSize.f90)
nca=50       #Number of Cells in Altgitude (no block division in alt) (if bui=1, will overwrite a line in the $hd copy of gitm/GITM2/src/ModSize.f90)
bui=0        #do you want this script to re-make GITM and re-make DART? Yes = 1, no = 0. Do 0, unless you changed some code or $est_f107 or $nbn, etc above.
f1t=150.0    #if you set $ts=1 above, what f107 would you like to use in that run (gets fed into clean.sh below). If you would like to use a real f107 (from f107.txt), uncomment lines 61 and 62 in clean.sh
f1m=150.0    #to do assimilation, you need to create an ensemble. The way it has been done so far is by creating several GITM instances with different f107's (sampled from normal distribution). What would you like the mean of f107 initial distribution to be? (gets fed to clean.sh below)
f1s=5.0      #what would you like the standard deviation of the initial distribution for f107 to be?
cfh=0.6      #horizontal cutoff, ie half-width of a Gaspari-Cohn polinomial for localizing observations, measured in horizontal great circle radians (not just latitude and longitude!). See section 8 of DART tutorial, slide 14 and cov_cutoff_mod.f90:line 142, and Gaspari and Cohn, 1999, QJRMS, eq 4.10 (will overwrite a line in the $rd copy of input.nml)
cfv=100000   #vert_normalization_height - like horizontal cutoff, but in vertical dimension, half-width of a Gaspari-Cohn polinomial. Measured in meters (will overwrite a line in the $rd copy of input.nml)
inf_v=1.01   #whole state inflation value - lambda on slide 6 of section 9 of DART tutorial and \lambda in eq (16) in gitm/GITM2/srcDoc/thermo.pdf and .tex and \gamma in Anderson and Anderson 1999 "A Monte Carlo Implementation of the Nonlinear Filtering Problem ...". This is a scaling up, so 1.01 corresponds to 1 percent increase over old values, 1.24 would be a 24% increase, etc. (will overwrite a line in the $rd copy of input.nml)
pr=1.0       #driver inflation value - constant variance of the driver as denoted by \sigma_i^2 in eq (18) in gitm/GITM2/srcDoc/thermo.pdf and .tex, measured in SFU^2 (will be communicated to GITM via sigma_i2.txt below, read into GITM in restart.f90)
aw=1800      #assimilation_period_seconds - ammount of simulated time between producing DART prior and posterior estimates, which requires GITM restarts, so making this period smaller will slow down the runs as writing GITM restarts takes time (will overwrite a line in the $rd copy of input.nml)
awt=300.0    #period between writing 3DALL/3DUSR files for the truth simulation (will overwrite line(s) in the $hd GITM2/srcData/UAM.in.test.DART.data_truth)
fk=1         #filter_kind in input.nml - see section 6 in DART tutorial - 1=EAKF, 2=EnKF, 3=Kernel, 4=particle, 5=rand samp, 6=kurtosis, 7=box kernel (will overwrite a line in the $rd copy of input.nml)

if [ $bui -eq 1 ]; then
    cd $hd/work
    cd ../GITM2
#change ModSize.f90 (GITM grid-size file)
    sed -i'.tmp' 's/  integer, parameter :: nLons =.*/  integer, parameter :: nLons ='$ncn'/' src/ModSize.f90
    sed -i'.tmp' 's/  integer, parameter :: nLats =.*/  integer, parameter :: nLats ='$nct'/' src/ModSize.f90
    sed -i'.tmp' 's/  integer, parameter :: nAlts =.*/  integer, parameter :: nAlts ='$nca'/' src/ModSize.f90
#change UAM.in (GITM input file) 
    sed -i'.tmp' 's/.*lons/'$nbn'    lons/' srcData/UAM.in.test.DART.data_ens
    sed -i'.tmp' 's/.*lats/'$nbt'    lats/' srcData/UAM.in.test.DART.data_ens
    sed -i'.tmp' 's/.*lons/'$nbn'    lons/' srcData/UAM.in.test.DART.data_truth
    sed -i'.tmp' 's/.*lats/'$nbt'    lats/' srcData/UAM.in.test.DART.data_truth
    sed -i'.tmp' 's/.*dt for output/'$awt'    dt for output/' srcData/UAM.in.test.DART.data_truth

    make || exit 1 


    cd ../work
    if [ $est_f107 -eq 0 ]; then #don't estimate f107
	sed -i'.tmp' 's/dist(k) =.*!changed.*/dist(k) = 1000.0 !changed by pbs_file script/' ../model_mod.f90 #every obs is very far away from f107 location (1000 radians is beyond physically meaningful 2pi radians on a sphere)
    else #estimate f107
	sed -i'.tmp' 's/dist(k) =.*!changed.*/dist(k) = 0 !changed by pbs_file script/' ../model_mod.f90 #every obs is exactly very close to the f107 location
    fi
    ./quickbuild.sh  || exit 2
fi

echo 'pbs_file.sh: start' `date` 

rm -rf $rd/${PBS_JOBNAME} #DELETES THE OLD DIRECTORY!
mkdir $rd/${PBS_JOBNAME} || exit 3 #makes a new directory on nobackup using the name of this job
cd $rd/${PBS_JOBNAME} 
cp -r $hd . || exit 4
cd gitm/work || exit 5

echo 'pbs_file.sh: past-cp' `date`

#change advance_model.csh to have the right number of processes
sed -i'.tmp' 's/@ nop = .* #changed.*/@ nop = '$nop' #changed by pbs_file script/' ../shell_scripts/advance_model.csh

#change input.nml (DART input file)
sed -i'.tmp' 's/ens_size.*=.*/ens_size = '$noe',/' input.nml   
sed -i'.tmp' 's/cutoff.*=.*/cutoff = '$cfh',/' input.nml
sed -i'.tmp' 's/vert_normalization_height.*=.*/vert_normalization_height = '$cfv',/' input.nml
sed -i'.tmp' 's/assimilation_period_seconds.*=.*/assimilation_period_seconds  = '$aw',/' input.nml
sed -i'.tmp' 's/filter_kind.*=.*/filter_kind = '$fk',/' input.nml
sed -i'.tmp' 's/inf_initial = 1.10, 1.0,/inf_initial = '$inf_v', 1.0,/' input.nml   

#choose the observations file
cp obs_seq.out.TEC obs_seq.out #TEC
#cp obs_seq.out.CHAMP obs_seq.out #CHAMP

pre=prespin_ #prefix for copies of prespinned folders

if [ $prs -eq 1 ]; then #make new prespin

    (./clean.sh $noe $nop $ts $f1t $f1m $f1s  >& output_c) || exit 5    
    for (( i = 1 ; i <= $[$noe + $ts]; i++ )) 
    do
	if [ $est_f107 -eq 0 ] ; then #don't estimate f107, so no need to inflate it, so keep the #DART=0 in all UAM.in's, nothing to change
	    cp -R advance_temp_e$i ${pre}advance_temp_e$i  #save a copy of the prespun directories
	else #estimate f107, so change #DART=0 flags in UAM.in's to appropriate settings:
	    
	    if [ $i -eq 1 ] ; then #first ens member is the master (1)
		sed -i'' 's/.*useDART (0=no, 1=master, 2=slave)/1  useDART (0=no, 1=master, 2=slave)/' advance_temp_e$i/UAM.in.truncated
		echo $noe > advance_temp_e$i/ens_size.txt
		echo $pr > advance_temp_e$i/sigma_i2.txt		
	    elif [ $i -eq $[$noe+1] ] ; then #truth is noe+1 (and no need for DART, hence 0)
		sed -i'' 's/.*useDART (0=no, 1=master, 2=slave)/0  useDART (0=no, 1=master, 2=slave)/' advance_temp_e$i/UAM.in.truncated   
	    else #other ens members are slaves (2)
		sed -i'' 's/.*useDART (0=no, 1=master, 2=slave)/2  useDART (0=no, 1=master, 2=slave)/' advance_temp_e$i/UAM.in.truncated
	    fi
	    
	    cp -R advance_temp_e$i ${pre}advance_temp_e$i #save a copy of the prespun directories

	fi
    done
    
else #prespin already exists
    
    ps=$rd/dart_test_ch/gitm/work #if you have already prespun the ensmble in one of the previous runs, tell me where it is
    
    cp ${ps}/perfect_ics . #only needed if you use perfect_model_obs
    cp ${ps}/filter_ics.* .
    for (( i = 1 ; i <= $[$noe + $ts]; i++ )) 
    do
	cp -r ${ps}/${pre}advance_temp_e$i advance_temp_e$i #copy the old prespin here and remove the prefix
	if [ $est_f107 -eq 0 ] ; then #don't estimate f107
	    echo ' ' #don't spread f107 means don't do anything as the prespin already had #DART=0
	else #estimate f107
	    if [ $i -eq 1 ] ; then #first ens member is the master (1)
		sed -i'' 's/.*useDART (0=no, 1=master, 2=slave)/1  useDART (0=no, 1=master, 2=slave)/' advance_temp_e$i/UAM.in.truncated
		echo $noe > advance_temp_e$i/ens_size.txt
		echo $pr > advance_temp_e$i/sigma_i2.txt		
	    elif [ $i -eq $[$noe+1] ] ; then #truth is noe+1 (and no need for DART, hence 0)
		sed -i'' 's/.*useDART (0=no, 1=master, 2=slave)/0  useDART (0=no, 1=master, 2=slave)/' advance_temp_e$i/UAM.in.truncated   
	    else #other ens members are slaves (2)
		sed -i'' 's/.*useDART (0=no, 1=master, 2=slave)/2  useDART (0=no, 1=master, 2=slave)/' advance_temp_e$i/UAM.in.truncated
	    fi  
	fi
    done    
fi

echo 'pbs_file.sh: past-ps' `date`

if [ $ts -eq 1 ]; then #if you requested a "truth spin"
#./run_p_fv.bash $noe $nop $vv || echo 1 POM fail  #the traditional DART way of doing it
    cd advance_temp_e$[$noe+1]
    cat ../../GITM2/srcData/UAM.in.test.DART.time_truth UAM.in > UAM.in.temporary #combine the data and time parts 
    mv UAM.in.temporary UAM.in
    echo " "             >> UAM.in 
    echo "#END"          >> UAM.in
    rm restartIN #remove symbolic link to restartOUT.out (where DART puts its restarts for GITM)
    rm UA/restartIN
    ln -s UA/restartOUT restartIN # and instead link it to the GITM's native restart location (restartOUT)
    ln -s restartOUT UA/restartIN

    tail -n $nop $PBS_NODEFILE > hf #create the hostfile for truth
    cat hf
    ( mpiexec -hostfile hf -n $nop GITM.exe < /dev/null > ens_log.txt ) & #run truth simulation on hf in background
    cd .. 
fi

echo 'pbs_file.sh: past-ts' `date`

(./run_f.bash $noe $[$noe*$nop] > output_f ) || exit 7 #nyx

echo 'pbs_file.sh: past-f' `date`

./obs_diag > obs_diag.out || exit 8

echo 'pbs_file.sh: past-od' `date`

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

