#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.
#
# WGL -- Though I am not great with shell scripts, I am hoping to update this
#          script for use with global WRF, which, most notably, does not 
#          boundary condition files.
#
# WGL -- should the following still be included, or is it peculiar to b.c.'s? ::
# If the ensemble mean assim_model_state_ic_mean is present in the CENTRALDIR,
# it is converted to a WRF netCDF format.

# Arguments are the process number of caller, the number of state copies
# belonging to that process, and the name of the filter_control_file for
# that process
set process = $1
set num_states = $2
set control_file = $3

set      myname = $0
set  CENTRALDIR = ..

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

set days_in_month = ( 31 28 31 30 31 30 31 31 30 31 30 31 )
   
set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'
set nonomatch
   
# give the filesystem time to collect itself
sleep 5

# Create a clean temporary directory and go there
${REMOVE} $temp_dir
mkdir -p  $temp_dir
cd        $temp_dir

# Each parallel task may need to advance more than one ensemble member.
# This control file has the actual ensemble number, the input filename,
# and the output filename for each advance.  Be prepared to loop and
# do the rest of the script more than once.
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3
while($state_copy <= $num_states)

   set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`

   set element = $ensemble_member
   
   # Shell script to run the WRF model from DART input.
   # where the model advance is executed as a separate process.
   
   echo "starting ${myname} for ens member $element at "`date`
   echo "CENTRALDIR is ${CENTRALDIR}"
   echo "temp_dir is ${temp_dir}"
   
   
   # Copy or link the required files to the temp directory
   
   ln -s ${CENTRALDIR}/input.nml .
   
   ln -s ${CENTRALDIR}/RRTM_DATA .
   ln -s ${CENTRALDIR}/LANDUSE.TBL .
   ln -s ${CENTRALDIR}/VEGPARM.TBL .
   ln -s ${CENTRALDIR}/SOILPARM.TBL .
   ln -s ${CENTRALDIR}/GENPARM.TBL .
   ln -s ${CENTRALDIR}/wrf.exe .
   ln -s ${CENTRALDIR}/gribmap.txt .

   # WGL:
   # Extra linked tables & data for global WRF (3.0) -- NOT all of these are
   #   always required, but it does no harm to link all of them so as to 
   #   avoid multiple working copies of advance_model.csh
   # These are always required 
   ln -s ${CENTRALDIR}/grib2map.tbl .
   # These are needed for use with module_ra_cam.F (I think)
   ln -s ${CENTRALDIR}/CAM_ABS_DATA .
   ln -s ${CENTRALDIR}/CAM_AEROPT_DATA .
   ln -s ${CENTRALDIR}/ozone.formatted .
   ln -s ${CENTRALDIR}/ozone_lat.formatted .
   ln -s ${CENTRALDIR}/ozone_plev.formatted .
   # These are needed for use with module_ra_gfdleta.F (I think)
   ln -s ${CENTRALDIR}/ETAMPNEW_DATA .
   ln -s ${CENTRALDIR}/tr49t67 .
   ln -s ${CENTRALDIR}/tr49t85 .
   ln -s ${CENTRALDIR}/tr67t85 .
   # These are needed for use with module_sf_urban.F (I think)
   # WRF versions before and including V3.0.1.1 use urban_param.tbl, 
   # V3.1 uses URBPARM.TBL - link to either that are found.)
   if ( -e ${CENTRALDIR}/urban_param.tbl ) then
      ln -s ${CENTRALDIR}/urban_param.tbl .
   endif
   if ( -e ${CENTRALDIR}/URBPARM.TBL ) then
      ln -s ${CENTRALDIR}/URBPARM.TBL .
   endif
   
   # nfile is required when using mpi to run wrf.exe
   # nfile is machine specific; ideally, it should be
   # constructed by the script advance_ens.csh
   
   hostname >! nfile
   hostname >> nfile
   ###ln -s  ${CENTRALDIR}/nfile$element nfile
   ${COPY} ${CENTRALDIR}/wrfinput_d0? .
                      # Provides auxilliary info not avail. from DART state vector
   
   if (  -e ${CENTRALDIR}/assim_model_state_ic_mean ) then
      ln -s ${CENTRALDIR}/assim_model_state_ic_mean dart_wrf_vector
      ${CENTRALDIR}/dart_to_wrf >& out.dart_to_wrf_mean
      ${COPY} wrfinput_d01 wrfinput_mean
   endif
   
   ${MOVE} ${CENTRALDIR}/$input_file dart_wrf_vector # ICs for run
   
   # Convert DART to wrfinput
   
   ${CENTRALDIR}/dart_to_wrf >& out.dart_to_wrf
   
   ${REMOVE} dart_wrf_vector
   
   # The program dart_to_wrf has created the file wrf.info.
   # Time information is extracted from wrf.info.
   
   set secday = `head -1 wrf.info`
   set targsecs = $secday[1]
   set targdays = $secday[2]
   set targkey = `echo "$targdays * 86400 + $targsecs" | bc`
   
   set secday = `head -2 wrf.info | tail -1`
   set wrfsecs = $secday[1]
   set wrfdays = $secday[2]
   set wrfkey = `echo "$wrfdays * 86400 + $wrfsecs" | bc`
   
   # If model blew up in the previous cycle, the member is now likely an outlier.
   
   if ( -e ${CENTRALDIR}/blown_${wrfdays}_${wrfsecs}.out ) then
      set MBLOWN = `cat ${CENTRALDIR}/blown_${wrfdays}_${wrfsecs}.out`
      set NBLOWN = `cat ${CENTRALDIR}/blown_${wrfdays}_${wrfsecs}.out | wc -l`
      set BLOWN = 0
      set imem = 1
      while ( $imem <= $NBLOWN )
         if ( $MBLOWN[$imem] == $element ) then
            set BLOWN = `expr $BLOWN \+ 1`
         endif
         set imem = `expr $imem \+ 1`
      end
   endif
   
   # START time from wrf.info file
   set cal_date    = `head -3 wrf.info | tail -1`
   set START_YEAR  = $cal_date[1]
   set START_MONTH = $cal_date[2]
   set START_DAY   = $cal_date[3]
   set START_HOUR  = $cal_date[4]
   set START_MIN   = $cal_date[5]
   set START_SEC   = $cal_date[6]
   
   # Initialize END time, but properly set it below
   set END_YEAR    = $cal_date[1]
   set END_MONTH   = $cal_date[2]
   set END_DAY     = $cal_date[3]
   set END_HOUR    = $cal_date[4]
   set END_MIN     = $cal_date[5]
   set END_SEC     = $cal_date[6]
   
   set MY_NUM_DOMAINS    = `head -4 wrf.info | tail -1`
   set ADV_MOD_COMMAND   = `head -5 wrf.info | tail -1`
   
   # Deal with Leap Year conventions
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
   
   ###################################################
   # Advance the model until target time is reached. #
   ###################################################
   
   # Integration length in seconds
   set INTERVAL_SS = `echo "$targkey - $wrfkey" | bc`

   # Set RUN time based on INTERVAL_SS -- do not use days, just hours and finer
   set RUN_DAYS    = 0
   set RUN_HOURS   = `expr $INTERVAL_SS \/ 3600`
   set REMAIN      = `expr $INTERVAL_SS \% 3600`
   set RUN_MINUTES = `expr $REMAIN \/ 60`
   set RUN_SECONDS = `expr $REMAIN \% 60`
   set INTERVAL_MIN = `expr $INTERVAL_SS \/ 60`
   
   # Set END time based on START time and RUN length
   set END_SEC = `expr $END_SEC \+ $INTERVAL_SS`
   while ( `expr $END_SEC - 60` >= 0 )
      set END_SEC = `expr $END_SEC \- 60`
      set END_MIN = `expr $END_MIN \+ 1`
      if ( `expr $END_MIN - 60` >= 0 ) then
         set END_MIN = `expr $END_MIN \- 60`
         set END_HOUR = `expr $END_HOUR \+ 1`
      endif
      if ( `expr $END_HOUR - 24` >= 0 ) then
         set END_HOUR = `expr $END_HOUR \- 24`
         set END_DAY  = `expr $END_DAY \+ 1`
      endif
      if ( `expr $END_DAY - $days_in_month[$END_MONTH]` > 0 ) then
         set END_DAY = 1
         set END_MONTH = `expr $END_MONTH \+ 1`
      endif
      if ( `expr $END_MONTH - 12` > 0 ) then
         set END_MONTH = 1
         set END_YEAR  = `expr $END_YEAR \+ 1`
   
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
   # Update time control entries in the WRF namelist.input:
   #-----------------------------------------------------------------------
   
   # Include RUN_DAYS too just in case!
   ${REMOVE} script.sed
   cat > script.sed << EOF
    /run_days/c\
    run_days                   = ${RUN_DAYS}
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
# set history interval equal to run interval to make sure you have
# a wrfoutput file at the end of the run interval
 /history_interval/c\
 history_interval           = ${INTERVAL_MIN}, ${INTERVAL_MIN}, ${INTERVAL_MIN}
#  dart_to_wrf is expecting only a single time per file
 /frames_per_outfile/c\
 frames_per_outfile         = 1, 1, 1,
EOF
   
   sed -f script.sed \
       ${CENTRALDIR}/namelist.input > namelist.input
   
   if ( -e rsl.out.integration ) then
      ${REMOVE} rsl.*
   endif
   
   # Set off WRF integration
   ${ADV_MOD_COMMAND} >>& rsl.out.integration
   ${COPY} rsl.out.integration ${CENTRALDIR}/wrf.out_${targdays}_${targsecs}_${element}  
   
   sleep 1
   
   set SUCCESS = `grep "wrf: SUCCESS COMPLETE WRF" rsl.* | cat | wc -l`
   if ($SUCCESS == 0) then
      if ($SUCCESS == 0) then
         echo $element >>! ${CENTRALDIR}/blown_${targdays}_${targsecs}.out
         echo "Model failure! Check file " ${CENTRALDIR}/blown_${targdays}_${targsecs}.out
         echo "for a list of failed elements, and check here for the individual output files:"
         echo " ${CENTRALDIR}/wrf.out_${targdays}_${targsecs}_elementnumber  "
         exit -1
      endif
   endif

   if ( -e ${CENTRALDIR}/extract ) then
      if ( $element == 1 ) then
         ls wrfout_d0${MY_NUM_DOMAINS}_* > wrfout.list
         if ( -e ${CENTRALDIR}/psfc.nc ) then
            ${COPY} ${CENTRALDIR}/psfc.nc .
         endif
         echo `cat wrfout.list | wc -l` | ${CENTRALDIR}/extract
         ${MOVE} psfc.nc ${CENTRALDIR}/.
      endif
   endif
   
   set dn = 1
   while ( $dn <= $MY_NUM_DOMAINS )
      ${MOVE} wrfout_d0${dn}_${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:${END_MIN}:${END_SEC} wrfinput_d0${dn}
      set dn = `expr $dn \+ 1`
   end
   
   ${REMOVE} wrfout*
   
   ##############################################
   # At this point, the target time is reached. #
   ##############################################
   
   # create new input to DART (taken from "wrfinput")
   ${CENTRALDIR}/wrf_to_dart >& out.wrf_to_dart
   
   ${MOVE} dart_wrf_vector ${CENTRALDIR}/$output_file

   set state_copy = `expr $state_copy \+ 1`
   set ensemble_member_line = `expr $ensemble_member_line \+ 3`
   set input_file_line = `expr $input_file_line \+ 3`
   set output_file_line = `expr $output_file_line \+ 3`
end



# Change back to working directory and get rid of temporary directory
cd ${CENTRALDIR}
echo ${REMOVE} ${temp_dir}
${REMOVE} ${temp_dir}

# Remove the filter_control file to signal completeion
# Is there a need for any sleeps to avoid trouble on completing moves here?
${REMOVE} $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

