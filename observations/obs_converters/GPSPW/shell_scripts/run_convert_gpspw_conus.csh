#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

set ddir = /glade/p/nmmm0024/syha/OBS_SEQ/GPSPW/data
set odir = /glade/p/nmmm0024/syha/OBS_SEQ/GPSPW
set wdir = `pwd`
set rdir = /glade/p/work/syha/DART/trunk/observations/GPSPW/work

cd $ddir
set fs = ( GPSPW_conus_201205*.nc )
#set fs = ( GPSPW_globe_20120528.nc GPSPW_globe_20120529.nc GPSPW_globe_2012053*.nc )
#           GPSPW_*_201206*.nc )
cd $wdir
if(! -e convert_gpspw) ln -s $rdir/convert_gpspw .
foreach f ( $fs )

  echo $f
  ln -s $ddir/$f gpspw_input.nc
  mv input.nml input.nml.temp
  set g = `echo $f | cut -d _ -f2`
  set t = `echo $f | cut -d _ -f3 | cut -d . -f1`
  #if($g == "conus") then
     set if_globe = .false.
  #else
  #   set if_globe = .true.
  #endif

  cat >! region.sed << EOF
 /global_data /c\
 global_data                 = ${if_globe}
EOF
 sed -f region.sed input.nml.temp > input.nml
 convert_gpspw > convert_gpspw.$g.$t.log
 
 set fout = obs_seq.gpspw.$g.$t
 if(! -e $fout) then
    echo Failed in convert_gpspw for $f.
    exit
 endif
 mv $fout $odir/
 ls -l $odir/$fout

end

exit 0


