#!/bin/bash  
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

#Aim: stitch together filter output to make diagnostic files

module load nco

echo -n "starting time "
date +"%T"

#Get new files, remove old ones
rm *.nc 
cp ../Posterior_Diag.nc .
cp ../Output/*.nc .

# Make copy the record dimension
ncpdq -O -a copy,time Posterior_Diag.nc Posterior_Diag.nc

export mean01='mean_copy_d01.nc'
export sd01='sd_copy_d01.nc'
export inf_mean01='post_inf_ic_new01_mean.nc'
export inf_sd01='post_inf_ic_new01_sd.nc'

ncecat -O -u copy $mean01 $mean01 &
ncecat -O -u copy $sd01 $sd01 &
ncecat -O -u copy $inf_mean01 $inf_mean01 &
ncecat -O -u copy $inf_sd01 $inf_sd01 &

wait

# rename variables
for var in U V MU PH T QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN H_DIABATIC U10 V10 T2 TH2 Q2 PSFC
do

   ncrename -v $var,$var'_d01' $mean01 &
   ncrename -v $var,$var'_d01' $sd01 &
   ncrename -v $var,$var'_d01' $inf_mean01 &
   ncrename -v $var,$var'_d01' $inf_sd01 &

   wait
 
   # concatenate the 4 files into Posterior_Diag
   # append
   ncrcat -A -v $var'_d01' $mean01 $sd01 $inf_mean01 $inf_sd01 Posterior_Diag.nc  
   echo -n "done variable " $var "time "
   date +"%T"   
done

# switch time and copy dimensions back
ncpdq -O -a time,copy Posterior_Diag.nc Posterior_Diag.nc
echo -n "finished " 
date +"%T"

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

