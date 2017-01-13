#!/bin/tcsh

set in = cesm.stdout.273578.locs
set out = cesm.stdout.273578.array

set obs = 600
set lines = 6

touch $out
set line = 1
set o = 1
while ($o <= $obs)
   set l = 1
   while ($l <= $lines)
      # perfect_model_obs: Processing observation        45  of    27000
      #{lon,lat} ob  =   9.331000000000E+01  0.000000000000E+00
      #{lon,lat}    12302 >   9.382917960675E+01 -4.363499758141E-15
      set list = `head -n $line $in | tail -n 1`	
      if ($l == 1) then
         set ob = $list[4]
      else if ($l == 2) then
         set list[3] = $ob
      endif
      if ($l > 1) then
         echo $list[3] $list[5-6] >> $out
      endif
      @ l++
      @ line++
   end
   echo "Finished ob $o"
   @ o++
end

exit
