#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

if ($#argv < 1) then
   echo "usage; check_cam num_ens_members"
   exit
endif 

set n = 1
while ($n <= $1)
     tail -40 cam_out_temp$n | grep 'END OF MODEL RUN' > /dev/null
     if ($status != 0) echo cam_out_temp$n finished abnormally
     @ n++
end

exit 0


