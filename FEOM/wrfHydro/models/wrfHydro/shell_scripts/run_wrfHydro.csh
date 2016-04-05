#!/bin/csh

set instance = $1
#set ppe = $1

# Touch wrfHydroStillWorking.dum is done at the top level of control to avoid lag
# resulting in too many cores being used.

if ( `readlink ../../wrf_hydro.exe | grep serial | wc -l` ) then 
    ../../wrf_hydro.exe >& wrfHydroOutput.${instance} 
else 
    mpirun -np 1 ../../wrf_hydro.exe >& wrfHydroOutput.${instance} 
endif 

\rm -f wrfHydroStillWorking.dum

exit 0
