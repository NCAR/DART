#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# for a model run. It assumes that there is ${PBS_O_WORKDIR}/WRF directory
# where boundary conditions files reside.

set days_in_month = ( 31 28 31 30 31 30 31 31 30 31 30 31 )

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Shell script to run the WRF model from DART input.

# set verbose

rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy or link the required files to the temp directory

mv ${PBS_O_WORKDIR}/assim_model_state_ic$element dart_wrf_vector # ICs for run
cp ${PBS_O_WORKDIR}/wrfinput_d0? . 
                   # Provides auxilliary info not avail. from DART state vector
ln -s ${PBS_O_WORKDIR}/input.nml .

ln -s  ${PBS_O_WORKDIR}/RRTM_DATA .
ln -s  ${PBS_O_WORKDIR}/LANDUSE.TBL .
ln -s  ${PBS_O_WORKDIR}/VEGPARM.TBL .
ln -s  ${PBS_O_WORKDIR}/SOILPARM.TBL .
ln -s  ${PBS_O_WORKDIR}/GENPARM.TBL .
ln -s  ${PBS_O_WORKDIR}/wrf.exe . 

hostname > nfile
hostname >> nfile
###ln -s  ${PBS_O_WORKDIR}/nfile$element nfile

# Convert DART to wrfinput

echo ".true." | ${PBS_O_WORKDIR}/dart_tf_wrf >& out.dart_to_wrf

set time = `head -1 wrf.info`
set targsecs = $time[1]
set targdays = $time[2]
set targkey = `expr $targdays \* 86400`
set targkey = `expr $targkey \+ $targsecs`

set time = `head -2 wrf.info | tail -1`
set wrfsecs = $time[1]
set wrfdays = $time[2]
set wrfkey = `expr $wrfdays \* 86400`
set wrfkey = `expr $wrfkey \+ $wrfsecs`

# Find all BC's file available and sort them with "keys".

#--1st, check if LBCs are "specified" (in which case wrfbdy files are req'd)
set SPEC_BC = `grep specified ${PBS_O_WORKDIR}/namelist.input | grep true | cat | wc -l`

if ($SPEC_BC > 0) then
   ls ${PBS_O_WORKDIR}/WRF/wrfbdy_*_$element > bdy.list
else
   echo ${PBS_O_WORKDIR}/WRF/wrfbdy_${targdays}_${targsecs}_$element > bdy.list
endif

echo ${PBS_O_WORKDIR}/WRF/wrfbdy_ > str.name
sed 's/\//\\\//g' < str.name > str.name2
set STRNAME = `cat str.name2`
set COMMAND = s/`echo ${STRNAME}`//

sed $COMMAND < bdy.list > bdy.list2
sed 's/_/ /g' < bdy.list2 > bdy.list
set num_files = `cat bdy.list | wc -l`
set items = `cat bdy.list`
set ifile = 1
set iday = 1
set isec = 2
while ( $ifile <= $num_files )
   set key = `expr $items[$iday] \* 86400`
   set key = `expr $key \+ $items[$isec]`
   echo $key >> keys
   @ ifile ++
   set iday = `expr $iday \+ 3`
   set isec = `expr $isec \+ 3`
end
set keys = `sort keys`

set date              = `head -3 wrf.info | tail -1`
set START_YEAR  = $date[1]
set START_MONTH = $date[2]
set START_DAY   = $date[3]
set START_HOUR  = $date[4]
set START_MIN   = $date[5]
set START_SEC   = $date[6]

set END_YEAR    = $date[1]
set END_MONTH   = $date[2]
set END_DAY     = $date[3]
set END_HOUR    = $date[4]
set END_MIN     = $date[5]
set END_SEC     = $date[6]

set MY_NUM_DOMAINS    = `head -4 wrf.info | tail -1`
set ADV_MOD_COMMAND   = `head -5 wrf.info | tail -1`

if ( `expr $END_YEAR \% 4` == 0 ) then
   set days_in_month[2] = 29
endif
if ( `expr $END_YEAR \% 100` == 0 ) then
   if ( `expr $END_YEAR \% 400` == 0 ) then
      set days_in_month[2] = 29
   else
      set days_in_month[2] = 28
   endif
endif

set ifile = 1
# Find the next BC's file available.

while ( $keys[${ifile}] <= $wrfkey )
   @ ifile ++
end

#################################################
# big loop here
#################################################

while ( $wrfkey < $targkey )

   set iday = `expr $keys[$ifile] \/ 86400`
   set isec = `expr $keys[$ifile] \% 86400`

# Copy the boundary condition file to the temp directory.
   cp ${PBS_O_WORKDIR}/WRF/wrfbdy_${iday}_${isec}_$element wrfbdy_d01

   if ( $targkey > $keys[$ifile] ) then
      set INTERVAL_SS = `expr $keys[$ifile] \- $wrfkey`
   else
      set INTERVAL_SS = `expr $targkey \- $wrfkey`
   endif
   set RUN_HOURS   = `expr $INTERVAL_SS \/ 3600`
   set REMAIN      = `expr $INTERVAL_SS \% 3600`
   set RUN_MINUTES = `expr $REMAIN \/ 60`
   set RUN_SECONDS = `expr $REMAIN \% 60`
   set INTERVAL_MIN = `expr $INTERVAL_SS \/ 60`

   @ END_SEC = $END_SEC + $INTERVAL_SS
   while ( $END_SEC >= 60 )
      @ END_SEC = $END_SEC - 60
      @ END_MIN ++
      if ($END_MIN >= 60 ) then
         @ END_MIN = $END_MIN - 60
         @ END_HOUR ++
      endif
      if ($END_HOUR >= 24 ) then
         @ END_HOUR = $END_HOUR - 24
         @ END_DAY ++
      endif
      if ($END_DAY > $days_in_month[$END_MONTH] ) then
         set END_DAY = 1
         @ END_MONTH ++
      endif
      if ($END_MONTH > 12 ) then
         set END_MONTH = 1
         @ END_YEAR ++
if ( `expr $END_YEAR \% 4` == 0 ) then
   set days_in_month[2] = 29
endif
if ( `expr $END_YEAR \% 100` == 0 ) then
   if ( `expr $END_YEAR \% 400` == 0 ) then
      set days_in_month[2] = 29
   else
      set days_in_month[2] = 28
   endif
endif
      endif
   end

   set END_SEC = `expr $END_SEC \+ 100`
   set END_SEC = `echo $END_SEC | cut -c2-3`
   set END_MIN = `expr $END_MIN \+ 100`
   set END_MIN = `echo $END_MIN | cut -c2-3`
   set END_HOUR = `expr $END_HOUR \+ 100`
   set END_HOUR = `echo $END_HOUR | cut -c2-3`
   set END_DAY = `expr $END_DAY \+ 100`
   set END_DAY = `echo $END_DAY | cut -c2-3`
   set END_MONTH = `expr $END_MONTH \+ 100`
   set END_MONTH = `echo $END_MONTH | cut -c2-3`

#-----------------------------------------------------------------------
# Create WRF namelist.input:
#-----------------------------------------------------------------------

rm -f script.sed
cat > script.sed << EOF
 /run_hours/c\
 run_hours                  = ${RUN_HOURS}
 /run_minutes/c\
 run_minutes                = ${RUN_MINUTES}
 /run_seconds/c\
 run_seconds                = ${RUN_SECONDS}
 /start_year/c\
 start_year                 = ${START_YEAR}, ${START_YEAR}, ${START_YEAR}
 /start_month/c\
 start_month                = ${START_MONTH}, ${START_MONTH}, ${START_MONTH}
 /start_day/c\
 start_day                  = ${START_DAY}, ${START_DAY}, ${START_DAY}
 /start_hour/c\
 start_hour                 = ${START_HOUR}, ${START_HOUR}, ${START_HOUR}
 /start_minute/c\
 start_minute               = ${START_MIN}, ${START_MIN}, ${START_MIN}
 /start_second/c\
 start_second               = ${START_SEC}, ${START_SEC}, ${START_SEC}
 /end_year/c\
 end_year                   = ${END_YEAR}, ${END_YEAR}, ${END_YEAR}
 /end_month/c\
 end_month                  = ${END_MONTH}, ${END_MONTH}, ${END_MONTH}
 /end_day/c\
 end_day                    = ${END_DAY}, ${END_DAY}, ${END_DAY}
 /end_hour/c\
 end_hour                   = ${END_HOUR}, ${END_HOUR}, ${END_HOUR}
 /end_minute/c\
 end_minute                 = ${END_MIN}, ${END_MIN}, ${END_MIN}
 /end_second/c\
 end_second                 = ${END_SEC}, ${END_SEC}, ${END_SEC}
 /history_interval/c\
 history_interval           = ${INTERVAL_MIN}, ${INTERVAL_MIN}, ${INTERVAL_MIN}
#  dart_tf_wrf is expecting only a single time per file
 /frames_per_outfile/c\
 frames_per_outfile         = 1, 1, 1,
EOF

 sed -f script.sed \
    ${PBS_O_WORKDIR}/namelist.input > namelist.input

# Update boundary conditions

   ${PBS_O_WORKDIR}/update_wrf_bc >& out.update_wrf_bc

   rm -f rsl.out.*

   ${ADV_MOD_COMMAND} >>& rsl.out.integration

   sleep 1

   set SUCCESS = `grep "wrf: SUCCESS COMPLETE WRF" rsl.* | cat | wc -l`
   if ($SUCCESS == 0) then
      echo $element >> $PBS_O_WORKDIR/blown_${targdays}_${targsecs}.out
   endif

   set dn = 1
   while ( $dn <= $MY_NUM_DOMAINS )
      mv -f wrfout_d0${dn}_${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:${END_MIN}:${END_SEC} wrfinput_d0${dn}
      @ dn ++
   end

   set START_YEAR  = $END_YEAR
   set START_MONTH = $END_MONTH
   set START_DAY   = $END_DAY
   set START_HOUR  = $END_HOUR
   set START_MIN   = $END_MIN
   set START_SEC   = $END_SEC
   set wrfkey = $keys[$ifile]
   @ ifile ++

end

#################################################
# end big loop here
#################################################

mv dart_wrf_vector dart_wrf_vector.input

# create new input to DART (taken from "wrfinput")
echo ".false." | ${PBS_O_WORKDIR}/dart_tf_wrf >& out.wrf_to_dart

mv -f dart_wrf_vector $PBS_O_WORKDIR/assim_model_state_ud$element

cd $PBS_O_WORKDIR
#rm -rf $temp_dir

exit
