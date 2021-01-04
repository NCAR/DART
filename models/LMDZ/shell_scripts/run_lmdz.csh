#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

limit stacksize unlimited
limit stacksize unlimited
limit datasize unlimited

source ../Control_File.csh

ln -sf $LMDZ_DEF_PATH/*.def .
ln -sf $LMDZ_DEF_PATH/$limit_file limit.nc
ln -sf $LMDZ_DEF_PATH/$gcm_exe .

rm used_*

../trans_time
set adv_date = `cat times | tail -1`
echo $adv_date
set hh = `echo $adv_date | cut -c12-13`
echo $hh
set ens_member = `cat element`
echo $ens_member

mv ../stok_paprs.dat_$ens_member  stok_paprs.dat

./$gcm_exe 

mv restart.nc start.nc
mv restartphy.nc startphy.nc

mv histhf.nc ../histhf_$ens_member.nc_$hh
mv histins.nc ../histins_$ens_member.nc_$hh

mv stok_paprs.dat ../stok_paprs.dat_$ens_member



exit 0


