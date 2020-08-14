#!/bin/csh
# 
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#>@todo FIXME ... this should use the advance_time program to manipulate the dates.
#>                advance_time gracefully handles leap years etc.

set month_days = (31 28 31 30 31 30 31 31 30 31 30 31)

set nmonth = 1
set imonth = 1

set year  = 2005
set m     = 4
while ( $imonth <= $nmonth )
   if ( $m > 12 ) then
       @ m = $m - 12
       @ year = $year + 1
   endif
  
   set iday = 2
   while ( $iday <= 2 )#$month_days[$m] )
      set month = `printf %02d $m`
      set day   = `printf %02d $iday`
      echo $year $month $day
      set obsdir = "$WORK/observations/syn/cesm2/ice-bridge/cice5_free_2005to2010/obs_seqs/aicen/err0.1/"

#     ls $obsdir/bootstrap/daily/obs_seqs/obs_seq.$year$month$day \
#        $obsdir/modis-tsfc/obs_seqs/obs_seq.$year-$month-$day-00000 > aice.tsfc.list

      ls $obsdir/obs_seq.aice?.${year}-${month}-${day}-00000  > cat.list

      sed "/filename_out/ c\      filename_out = '/$obsdir/obs_seq.$year-$month-$day-00000'" input.nml>temp

      mv temp input.nml
      ./obs_sequence_tool
      @ iday = $iday + 1
   end
   @ imonth = $imonth + 1
   @ m      = $m + 1
end

exit 0

