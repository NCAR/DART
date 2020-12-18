#!/bin/csh
# 
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#>@todo FIXME ... this should use the advance_time program to manipulate the dates.
#>                advance_time gracefully handles leap years etc.

set month_days = (31 28 31 30 31 30 31 31 30 31 30 31)
set case  = ( osse_inf_loc0.01_sst osse_inf_loc0.01_pw_sst osse_inf_loc0.01_tw_sst )

set nmonth = 36
set imonth = 1

set year  = 2001
set m     = 1
while ( $imonth <= $nmonth )
   if ( $m > 12 ) then
       @ m = $m - 12
       @ year = $year + 1
   endif
  
   set iday = 1
   while ( $iday <= $month_days[$m] )
      set month = `printf %02d $m`
      set day   = `printf %02d $iday`
      echo $year $month $day

      ls $SCRATCH/$case[1]/Obs_seqs/cice.obs_seq.$year-$month-$day-00000.final \
         $SCRATCH/$case[2]/Obs_seqs/cice.obs_seq.$year-$month-$day-00000.final \
         $SCRATCH/$case[3]/Obs_seqs/cice.obs_seq.$year-$month-$day-00000.final  > list1.txt

      ./obs_common_subset
      @ iday = $iday + 1
   end
   @ imonth = $imonth + 1
   @ m      = $m + 1
end

exit 0

