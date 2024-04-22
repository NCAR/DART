#!/bin/bash

# Directories of the truth run
# and the routelink file here:
truth_d=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/
route_l=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/RouteLink.nc
gageind=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/gage_ind_small.txt
gageids=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/gauge_ids.txt
cd $truth_d

echo ' '

rm combo_IDs || true

ncdump -f F -v gages RouteLink.nc | grep "  01" | cut -d, -f1 | tr '"' ' ' | sed 's/^[ \t]*//;s/[ \t]*$//' > gauges_rl
ncdump -f F -v gages RouteLink.nc | grep "  01" | cut -d, -f3 | tr ")" " " > indxes_rl

num_gauges_rl=`cat gauges_rl | wc -l`
echo $num_gauges_rl

k=1
while read -r line; do
    gauges[$k]=$line
    let 'k+=1'
done < "gauges_rl"

k=1
while read -r line; do
    indxes[$k]=$line
    let 'k+=1'
done < "indxes_rl"

# All gauges available in the RL file
for k in `seq 1 $num_gauges_rl`; do
    echo $k,${gauges[$k]},${indxes[$k]} >> combo_IDs
    echo Gauge: $k, ID: ${gauges[$k]}, Feature: ${indxes[$k]} 
done


# Fetch Gauges/Links/
# If both link file and gauge file are given, priority 
# is for the gauge IDs: 
if   [[ -f "${gageids}" ]]; then 
     # User provided set of gauges 
     while read -r line; do
         grep $line combo_IDs | cut -d, -f3 >> tmp        
     done < "${gageids}"
     echo -e "\nThe truth will be computed at user-provided gauge ID locations."

elif [[ -f "${gageind}" ]]; then
     # User provided set of links
     cp ${gageind} ${truth_d}/tmp
     echo -e "\nThe truth will be computed at user-provided index locations."

else
     # Couldn't find files that contain either gauge IDs or links 
     mv indxes_rl tmp
     echo -e "\nThe gage locations from the RouteLink file are used to form the truth."
fi
