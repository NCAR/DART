#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
From chriss@ucar.edu Fri Dec 20 15:10:47 2002
Date: Thu, 19 Dec 2002 16:18:10 -0700 (MST)
From: Chris Snyder <chriss@ucar.edu>
To: Jeffrey Anderson <jla@ucar.edu>, Dale Barker <dmbarker@ucar.edu>
Subject: Bill's enkf script

#
#  script for execution of ensemble Enkf filter.
#  members are integrated on various remote machines.

#--- Variables ------
#set verbose

echo ' begin enkf run '
set rundir = /data5/chriss/Skam/Test/DA_run/Run
                                # in ./Truth and ./Member.*** dirs
set datadir =  $rundir/Archive  # Data copied to here
set filterdir = $rundir/Filter  # Filter runs in here
set sounding = sq10

set t0     = 0      # Ensemble initialized at this time
set t1     = 50     # First obs available at this time
set deltat = 50     # Interval between obs (or could specify string of times)
set Ncyc   = 18     # Number of assimilation cycles
set Nparallel = 14  # Number of members to run in parallel

set Nens   = 80     # Number of ensemble members
set programs = ( cloud3d cloud3d_enkf init \
                 init_perfect )
set datafiles = ( 3d.run 3d.params \
                  fip.make_members fip.apply_filter \
                  sq10.sound cloud3d_extract_obs.ncl \
                  cloud3d_input_truth \
                  wrf_user_fortran_enkf.so wrf_user_cloud_lf.so )

set machines = ( service01 service02 service03 service04 service05  \
                 service06 service07 service08 service09 service10  \
                 service11 service12 service13 service14 )

set master   = paloverde

# UNITS  -> number of time steps

rm -f forecast_abort

# list of assimilation times

echo ' set times '

set times = `pad $t1`
set t = `pad $t1`
set iloop = 1 

@ max = $Ncyc - 1 

while ( $iloop <= $max )  
  @ t = $t + $deltat
  set t = `pad $t`
  set times = ( $times $t ) ; @ iloop++
end

echo ' times are ' $times

# list of members
set ensmems = () 
set iloop = 1 
while ( $iloop <= $Nens )  
  set iens = `pad $iloop`
  set ensmems = ( $ensmems $iens ) ; @ iloop++
end
set ensm = $Nens
@ ensm++
set ensmean = `pad $ensm`

echo ' list of members ' $ensmems
echo ' mean is         ' $ensmean

#--- Set up -----------

cd $rundir

#  directories for the parallel (remote) integration of members

mkdir Need_to_run
mkdir Running
mkdir Finished

foreach prog ($programs[*])
  if (! -x $prog) then
    echo ' missing program ' $prog
    exit
  endif
end

foreach file ($datafiles[*])
  if (! -e $file) then
    echo ' missing data file ' $file
    exit
  endif
end

if(! -d $datadir) then
  mkdir $datadir
endif

if(! -d $filterdir) then
  mkdir $filterdir
endif

if(! -d Truth) then
  mkdir Truth
endif

foreach iens ($ensmems[*])
  mkdir Member.$iens  # Make dirs for each member
  cp 3d.run Member.$iens/3d.run
  cp cloud3d Member.$iens/cloud3d
end

 echo " initializing truth, members "

 cp $sounding.sound Truth/$sounding.sound

 cp 3d.params Truth/3d.params
 cp 3d.run Truth/3d.run

#  this makes the base for the ensemble members
   cd Truth; $rundir/init > &out.init

   mv $sounding.r00000.g01 $filterdir/cloud3d_enkf_input_0001

#  this makes the initial conditions for the truth run
#   $rundir/init_perfect.exe > &out.init_perfect
   cp $rundir/cloud3d_input_truth $sounding.r00000.g01

 cd $rundir

#--- set initial ensemble members
#

#  ********* 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  *********  this member initialization uses the cloud3d
#  *********  initialization code to generate the members
#

foreach iens ($ensmems[*])
  init > out.init
  mv $sounding.r00000.g01 Member.$iens/$sounding.r00000.g01
end

#  ********* 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#  *******  this member initialization puts 
#  *******  a random perturbation on a single state
#
# cp  fip.make_members $filterdir/filter_input_parameters
# cd $filterdir
#   $rundir/cloud3d_enkf > &out.cloud_enkf_init
#
#   foreach iens ($ensmems)
#     cp cloud3d_enkf_output_$iens $datadir/a$t0.c000.m$iens  # archive IC for each member
#                                                             # convention: a = analysis
#                                                             # time (no padding???)
#                                                             # cycle #  (3 digits, padded)
#                                                             # member # (4 digits, padded)
#     mv cloud3d_enkf_output_$iens $rundir/Member.$iens/$sounding.r00000.g01
#
#   end
#
  cd $rundir
 
cp fip.apply_filter  $filterdir/filter_input_parameters # copy filter input
cp cloud3d_extract_obs.ncl     Truth/
cp wrf_user_fortran_enkf.so     Truth/
cp wrf_user_cloud_lf.so     Truth/

mkdir Machine_outputs

#--- Cycle EnKF -------

set told = t0; set icycle = 1;
foreach t ($times)

  set c = `pad $icycle 3` # cycle number for archived files

  echo "***************************************************"
  echo "      beginning filter cycle " $c " at " `date`
  echo "***************************************************"

  #--- forecasts

  # set up the Need_to_run directory

    foreach iens ($ensmems[*])
      touch Need_to_run/Member.$iens
    end
    rm -f Running/*
    rm -f Finished/*

   #  for each available machine, copy over the necessary script and run it

  foreach machine ($machines[1-$Nparallel])
#    rcp Machine.scripts/$machine.script $machine":"$machine.script
    rcp Machine.scripts/service01.script $machine":"$machine.script
#    if(! -e pad ) rcp Machine.scripts/pad $machine":"pad

#    echo ' starting ensembles on machine ' \
#          $machine ' for cycle ' $c >>  out.$machine
    (rsh $machine $machine.script \
     $machine $master $rundir $Nens cloud3d 3d.run sq10.r00000.g01 \
             >> & Machine_outputs/out.$machine.c$c) &
  end

  #--- integrate perfect model, get new obs, we'll do this locally for now

    cd Truth

    echo " starting truth, cycle " $c " at " `date`
    $rundir/cloud3d > & cloud3d_output_$c

      if ($status != 0) then
        echo " model abort in truth run, error exit for enkf "
        exit
      endif

    mv $sounding.r00000.g01 cloud3d_input.$c
    mv $sounding.r* cloud3d_output  # move model output 
                                    # for input to obs program
    echo " finished truth, cycle " $c  " at " `date`
    echo " extracting obs, cycle " $c 
    ncl < cloud3d_extract_obs.ncl > out_obs_$c
    echo " finished obs, cycle   " $c  " at " `date`
    mv wrf_obs $filterdir/cloud3d_enkf_obs
    mv cloud3d_output $sounding.r00000.g01

    cd $rundir

  wait #  wait for the members to be integrated

#  push the forecast files to the Filter, and archive forecasts

  foreach iens ($ensmems)
    cd Member.$iens

    # check for model blowup
    if( -e abort.member.$iens ) then
      echo " member " $iens " integration aborted "
      echo " error exit "
      exit
    endif

    \rm -f $sounding.r00000.g01  # delete start data for member,
                                 # its already been archived
    cat Member.$iens.timing out.member.$iens > output.m$iens.c$c
    rm -f Member.$iens.timing out.member.$iens
    cp $sounding.r*.g01 $datadir/f$t.c$c.m$iens    # archive forecasts
    mv $sounding.r*.g01 $filterdir/cloud3d_enkf_input_$iens 
                                             # placed member in filter dir

    cd $rundir
  end


  #--- EnKF

  cd $filterdir

  echo " starting filter for cycle " $c  " at " `date`
  $rundir/cloud3d_enkf > &out.cloud3d_enkf.c$c

      if ($status != 0) then
        echo " filter abort, error exit for enkf "
        exit
      endif

  echo " finished filter for cycle " $c  " at " `date`

  foreach iens ($ensmems)
# don't archive analysis, we can recreate from forecasts
#    cp cloud3d_enkf_output_$iens $datadir/a$t.c$c.m$iens  # archive analyses
       # analysis as IC for next f/c
    mv  cloud3d_enkf_output_$iens $rundir/Member.$iens/$sounding.r00000.g01 
  end

  mv cloud3d_enkf_output_$ensmean $datadir/a$t.c$c.m$ensmean
  mv fort.10 mean_output.$c
  \rm -f cloud3d_enkf_input_*

  cd $rundir

  set told = t  # save time for use in next cycle
  @ icycle++    # increment cycle number

  echo "***************************************************"
  echo "      finished filter cycle " $c " at " `date`
  echo "***************************************************"

end

  echo "***************************************************"
  echo "      Enkf completed " $c " at " `date`
  echo "***************************************************"
