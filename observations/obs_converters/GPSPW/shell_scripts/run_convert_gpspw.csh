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
set fs = ( GPSPW_globe_201206*.nc )
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
  #   set if_globe = .false.
  #else
     set if_globe = .true.
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
 set gday = `grep OBS_DATE convert_gpspw.$g.$t.log | awk '{print $2}'`
 mv input.nml input.nml.temp

 foreach ih ( 00 06 12 18 )
 set fnew = $odir/${fout}${ih}
 set   is = `expr $ih \* 3600`
 cat >! date.sed << EOF
  /filename_seq /c\
   filename_seq         = '$fout',
  /filename_out /c\
   filename_out         = '$fnew',
  /first_obs_days /c\
   first_obs_days       = $gday,
  /first_obs_seconds /c\
   first_obs_seconds    = $is,
  /last_obs_days /c\
   last_obs_days        = $gday,
  /last_obs_seconds /c\
   last_obs_seconds     = $is,
EOF
 sed -f date.sed input.nml.temp > input.nml
 obs_sequence_tool > obs_seq_tool.$g.$t$ih.log 
 if(! -e $fnew) then
    echo Failed in obs_sequence_tool for $fnew.
    exit
 endif
 ls -l $fnew
 end

end

exit 0


