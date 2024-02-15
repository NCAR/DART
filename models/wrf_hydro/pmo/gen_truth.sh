#!/bin/bash

# Directories of the truth run, gauges indices, gauge IDs
# and the routelink file here:
ref_dir=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/
route_l=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/RouteLink.nc
gageind=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/gage_ind.txt
gageids=/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/gage_ids.txt

cd ${ref_dir}

echo ' '

# time stamps: OSSE start to OSSE end
osse_s=2019-06-01_00
osse_e=2019-06-01_23

# hydro restart file 
hydroL=HYDRO_RST.
hydroR=_DOMAIN1

# Truth file name
truth_x=pmo_drb_truth.nc
truth_d=pmo_drb_truth_all_gages.nc
truth_g=pmo_drb_truth_des_gages.nc

rm -f $truth_x $truth_d $truth_g

# Form a list of all needed files
ls -1 ${hydroL}*${hydroR} > majorlist.txt

line1=`grep -n $osse_s majorlist.txt | head -n 1 | cut -d: -f1`
line2=`grep -n $osse_e majorlist.txt | head -n 1 | cut -d: -f1`

sed -n "${line1},${line2}p" majorlist.txt > newlist.txt

nhoura=`cat newlist.txt | wc -l`
nhourb=`echo "($nhoura - 1)" | bc -l`
timesq=`seq 0 $nhourb`

# Remove all variables and keep qlink1
f=0

for file in `cat newlist.txt`; do

    let "f+=1"

    echo $file, Cycle: $f

    ex=`printf '%05.f' $f`

    # Extract qlink1 only 
    ncks -O -v qlink1 $file member${ex}.nc
   
    # Rename variable and dimension
    ncrename -O -d links,feature_id -v qlink1,streamflow member${ex}.nc member${ex}.nc

    # Add record time dimension
    ncecat -O -u time member${ex}.nc member${ex}.nc
done

# Concatenate all files
ncrcat -F -O -h -H member?????.nc $truth_x

rm member*.nc || true

ncap2 -O -s 'streamflow@units="m3 s-1";streamflow@long_name="River Flow";streamflow@grid_mapping="crs"; streamflow@missing_value=-9999.f' $truth_x $truth_x
ncatted -O -a cell_methods,streamflow,d,, $truth_x $truth_x

# Add time variable
ncap2 -O -s 'time=array(0,1,$time)' $truth_x $truth_x
ncap2 -O -s 'time@long_name="valid output time";time@standard_name="time";time@units="hours since ${osse_s}";time@calendar="proleptic_gregorian"' $truth_x $truth_x
ncatted -O -a units,time,o,c,"hours since ${osse_s}:00" $truth_x

# Bring in the feature IDs from RL
ncks -A -v link $route_l $truth_x
ncrename -O -v link,feature_ids $truth_x $truth_x  

# Clean-up some global attributes
ncatted -O -a Restart_Time,global,d,, $truth_x
ncatted -O -a Since_Date,global,d,, $truth_x
ncatted -O -a his_out_counts,global,d,, $truth_x
ncatted -O -a featureType,global,a,c,"timeSeries" $truth_x
ncatted -O -a station_dimension,global,a,c,"feature_id" $truth_x
ncatted -O -a NCO,global,d,, $truth_x
ncatted -O -a history,global,d,, $truth_x

echo -e "\n** Created reference truth trajectory: ${ref_dir}${truth_x}\n"


# Still need to subset based on gages
ncdump -f F -v gages RouteLink.nc | grep "  01" | cut -d, -f1 | tr '"' ' ' | sed 's/^[ \t]*//;s/[ \t]*$//' > gauges_rl
ncdump -f F -v gages RouteLink.nc | grep "  01" | cut -d, -f3 | tr ")" " " > indxes_rl

num_gauges_rl=`cat gauges_rl | wc -l`
echo -e "** The number of gauges in the Route Link file is: $num_gauges_rl\n"

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

sleep 0.5

# All gauges available in the RL file
for k in `seq 1 $num_gauges_rl`; do
    echo $k,${gauges[$k]},${indxes[$k]} >> combo_IDs
    echo Gauge: $k, ID: ${gauges[$k]}, Index: ${indxes[$k]} 
done

# Fetch Gauges/Links/
# If both link file and gauge file are given, priority 
# is for the gauge IDs: 
if   [[ -f "${gageids}" ]]; then
     # User provided set of gauges 
     while read -r line; do
         grep $line combo_IDs | cut -d, -f3 >> tmp
     done < "${gageids}"
    
     num_gauges_des=`cat ${gageids} | wc -l`

     echo -e "\n** The truth will be computed at user-provided gauge ID locations only."

elif [[ -f "${gageind}" ]]; then
     # User provided set of links
     cp ${gageind} ${ref_dir}/tmp

     num_gauges_des=`cat ${gageind} | wc -l`

     echo -e "\n** The truth will be computed at user-provided index locations only."

else
     # Couldn't find files that contain either gauge IDs or links 
     cp indxes_rl tmp
     
     num_gauges_des=${num_gauges_rl}

     echo -e "\n** The gage locations from the RouteLink file are used to form the truth."
fi

# Permutate the record to feature_id instead of time
ncpdq -O -a feature_id,time $truth_x new_pmo.nc

echo -e "\n** Creating truth file at the desired gauges: ${ref_dir}$truth_g"

# Create individual files at the gage locations
cc=0
for i in `cat tmp`; do
   let 'cc+=1'
   
   fid=$(($i-1))
   ncks -O -v feature_ids,time,streamflow -d feature_id,$fid new_pmo.nc test_${cc}.nc
done

# Now, concatenate the resulting files
ncrcat -O test_*.nc $truth_g 
ncpdq -O -a time,feature_id $truth_g $truth_g
ncatted -O -a history,global,d,, $truth_g

# Check if we need to provide the truth for all gauges
if [[ $num_gauges_rl != $num_gauges_des ]]; then
   cp indxes_rl tmp
 
   cc=0
   for i in `cat tmp`; do
      let 'cc+=1'
   
      fid=$(($i-1))
      ncks -O -v feature_ids,time,streamflow -d feature_id,$fid new_pmo.nc test_a${cc}.nc
   done

   ncrcat -O test_a*.nc $truth_d
   ncpdq -O -a time,feature_id $truth_d $truth_d
   ncatted -O -a history,global,d,, $truth_d
else
   cp $truth_g $truth_d
fi

rm test*.nc new_pmo.nc tmp combo_IDs indxes_rl gauges_rl || true

echo -e "\n          ##### Done #####"
