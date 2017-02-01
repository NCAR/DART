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
cp ../Prior_Diag.nc .
cp ../Output/*.nc .

# Make copy the record dimension
ncpdq -O -a copy,time Prior_Diag.nc Prior_Diag.nc

export mean='prior_diag_mean01.nc'
export sd='prior_diag_sd01.nc'
export inf_mean='prior_diag_inf_mean01.nc'
export inf_sd='prior_diag_inf_sd01.nc'

# create time dimension in POP output
ncecat -O -u time $mean $mean &
ncecat -O -u time $sd $sd &
ncecat -O -u time $inf_mean $inf_mean &
ncecat -O -u time $inf_sd $inf_sd &

wait

# create copy record dimension in POP output
ncecat -O -u copy $mean $mean &
ncecat -O -u copy $sd $sd &
ncecat -O -u copy $inf_mean $inf_mean &
ncecat -O -u copy $inf_sd $inf_sd &

wait

## rename variables
for var in 'SALT' 'TEMP' 'UVEL' 'VVEL' 'PSURF' 
do

   ncrename -v $var'_CUR',$var $mean &
   ncrename -v $var'_CUR',$var $sd &
   ncrename -v $var'_CUR',$var $inf_mean &
   ncrename -v $var'_CUR',$var $inf_sd &

   wait
 
   # concatenate the 4 files into Prior_Diag
   # append
   ncrcat -A -v $var $mean $sd $inf_mean $inf_sd Prior_Diag.nc  
   echo -n "done variable " $var "time "
   date +"%T"   
done

# switch time and copy dimensions back
ncpdq -O -a time,copy Prior_Diag.nc Prior_Diag.nc
echo -n "finished " 
date +"%T"

#mv Prior_Diag.nc Prior_Diag.nc.full

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

