#!/bin/csh -f
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

set argv = (`getopt Hbehmsf:t: $*`)

#-----------------------------------------------------------------------

 set sep =  
 set format = standard

 set rec = tail
 set help = 0

 set hours   = 0
 set minutes = 0
 set seconds = 0

#-----------------------------------------------------------------------

while ("$argv[1]" != "--")
    switch ($argv[1])
        case -H:
            set help = 1; breaksw
        case -b:
            set rec = head; breaksw
        case -e:
            set rec = tail; breaksw
        case -h:
            set hours = 1; breaksw
        case -m:
            set hours = 1; set minutes = 1; breaksw
        case -s:
            set hours = 1; set minutes = 1; set seconds = 1; breaksw
        case -f:
            set format = $argv[2]; shift argv; breaksw
        case -t:
            set sep = $argv[2]; shift argv; breaksw
    endsw
    shift argv
end
shift argv

#  --- help output ---

 if ( $help ) then
   cat << END
   time_stamp.csh [ -behms -f format -t separator ]

      -H help (no execution)
      -b beginning date
      -e ending date (default)
      -h hours
      -m hours,minutes
      -s hours,minutes,seconds
      -f format=standard(default),european,digital
      -t separator (default=blank)

END
   exit
 endif

#  --- check format ---

 if ( $format != "standard" &&  \
      $format != "european" && $format != "digital" ) then
    echo ERROR invalid format
    exit (4)
 endif
 
 set hsep = $sep
 if ( $format == "standard" || $format == "european" ) set hsep = h

#-----------------------------------------------------------------------

 if ( -e time_stamp.out ) then
     set time_stamp = `$rec -1 time_stamp.out`

     set month_name = `echo $time_stamp[7] | tr "[A-Z]" "[a-z]"`
     set  month_num = $time_stamp[2]
     if ($month_num < 10) set month_num = "0"$time_stamp[2]

#    ---- day can have more than 2 digits ----

     set digits = 2
     if ( $month_name == "day" ) set digits = 4

     set   day = $time_stamp[3]
     set   day_num = $day
     set i = 1
     set x = 1
     while ( $i < $digits )
        @ i++
        @ x = $x * 10
        if ( $day < $x ) set day_num = "0"$day_num
     end

#    ---- hours,min,sec can have only 2 digits ----

     set  hour_num = $time_stamp[4]
     if ($hour_num < 10) set hour_num = "0"$time_stamp[4]

     set  min_num = $time_stamp[5]
     if ($min_num < 10) set min_num = "0"$time_stamp[5]

     set  sec_num = $time_stamp[6]
     if ($sec_num < 10) set sec_num = "0"$time_stamp[6]

#    ---- create date label ----

     set date_name

     if ( $format == "standard" ) then
        if ( $month_name != "day" ) set date_name = $time_stamp[1]
        set date_name = $date_name$month_name$day_num
     else if ( $format == "european" ) then
        set date_name = $day_num$month_name$time_stamp[1]
     else if ( $format == "digital" ) then
        set date_name = $time_stamp[1]$sep$month_num$sep$day_num
     endif

        if ( $hours   ) set date_name = $date_name$hsep$hour_num
        if ( $minutes ) set date_name = $date_name$sep$min_num
        if ( $seconds ) set date_name = $date_name$sep$sec_num

 else
#    --- dummy values ---
     set month_name = "xxx"
     set date_name  = "no_time_stamp"
 endif

     echo $date_name

