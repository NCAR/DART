#!/bin/tcsh

# Tally the number of obs at each unique time.  
# Reads a file created from
# $ set pre_day  = the day before the one we want to summarize
# $ set post_day = the day after the one we want to summarize
# $ grep -B 2 OBS obs_seq.out | \
#   grep -v -e OBS -e '^[ ]*-' -e '--' -e first -e $pre_day \
#        -e $post_day -e COSMIC  >! obs_seq.out.times


#Output looks like
# sec   day    # obs
# 86336 153035 131
# 10001 153036 126

set ofile = $1

set line = `wc -l $ofile`
set nlines = $line[1]

if (-f $ofile.uniq) rm $ofile.uniq
touch $ofile.uniq

# KR.times
#  86318     153035
#  86336     153035
#  10001     153036
#  10007     153036
# KR.times.uniq
# 86336 153035 131
# 10001 153036 126
# 10001 153036 0

set l = 1
set line = ( `sed -n ${l}p $ofile` )
set sec = $line[1]
set day = $line[2]
set obs = 1
# echo $sec $day $obs
@ l++

while ($l < $nlines)
# while ($l < 500)
   set line = ( `sed -n ${l}p $ofile` )
   if ($line[1] == $sec  &&  $line[2] == $day) then
      @ obs++
   else
# The original output consisted of the line just read and the previous $obs.
      echo $sec $day $obs >> $ofile.uniq
      set sec = $line[1]
      set day = $line[2]
      set obs = 0
#       echo "   $sec $day $obs"
   endif
   @ l++
end
