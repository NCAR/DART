#!/bin/tcsh

# script to extract the platforms (fake obs err var) of a list of times.
# The list of times could come from uniq_times.csh.

set tfile = $1
set ofile = $2
set day = $3

set line = `wc -l $tfile`
set ntimes = $line[1]

if (-f $tfile.platforms) mv $tfile.platforms $tfile.platforms.$$
touch $tfile.platforms

set t = 1
while ($t <= $ntimes)
# while ($t <= 10)
   set line = ( `sed -n ${t}p $tfile` )
   set sec = $line[1]

#    set line = `grep -m 1 -A 1 "^ $sec" $ofile | tail -n 1` 
   set line = `grep -m 1 -A 1 "^[ ]*$sec     $day" $ofile | tail -n 1` 
   echo "time = $sec, platform = $line" >> $tfile.platforms

   @ t++
end
