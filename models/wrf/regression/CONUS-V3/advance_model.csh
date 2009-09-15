#!/bin/csh
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

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# for a model run. It assumes that there is ${CENTRALDIR}/WRF directory
# where boundary conditions files reside.

# If the ensemble mean assim_model_state_ic_mean is present in the CENTRALDIR,
# it is converted to a WRF netCDF format.
# It is then used in update_wrf_bc the calculate the deviation from the mean.
# This deviation from the mean is then added at the end of the interval to
# calculate new boundary tendencies. The magnitude of the perturbation added
# at the end of the interval is controlled by infl. The purpose is to increase
# time correlation at the lateral boundaries.

# Arguments are the process number of caller, the number of state copies
# belonging to that process, and the name of the filter_control_file for
# that process
set process = $1
set num_states = $2
set control_file = $3
# Setting to 1 saves output files from the ensemble mean only, while setting to
# larger numbers will save all member output files <= to this value
set save_elements = 1

set      myname = $0
set  CENTRALDIR = /ptmp/romine/work/cv3work
set  WRFOUTDIR  = /ptmp/romine/work/cv3work/wrfout

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

set days_in_month = ( 31 28 31 30 31 30 31 31 30 31 30 31 )
   
set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'
   
# give the filesystem time to collect itself
sleep 5

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

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
   
   set infl = 0.0
   
   # Shell script to run the WRF model from DART input.
   # where the model advance is executed as a separate process.
   
   # echo "starting ${myname} for ens member $element at "`date`
   # echo "CENTRALDIR is ${CENTRALDIR}"
   # echo "temp_dir is ${temp_dir}"
   
   
   # Copy or link the required files to the temp directory
   
   ln -sf ${CENTRALDIR}/input.nml .
   ln -sf ${CENTRALDIR}/wrf.exe   .
 
   # These are all WRF aux files; one by one, link them into the
   # temp running directory.  It is currently NOT an error if the
   # file does not exist -- this makes the list compatible with different
   # versions of WRF (e.g. v2 vs v3) but be aware that wrf may fail to run
   # if the required files are not found.

#   set wrf_aux_file_list = ' \
#      CAM_ABS_DATA \
#      CAM_AEROPT_DATA \
#      ETAMPNEW_DATA \
#      ETAMPNEW_DATA_DBL \
#      GENPARM.TBL \
#      grib2map.tbl \
#      gribmap.txt \
#      LANDUSE.TBL \
#      ozone.formatted \
#      ozone_lat.formatted \
#      ozone_plev.formatted \
#      RRTM_DATA \
#      RRTM_DATA_DBL \
#      RRTMG_LW_DATA \
#      RRTMG_LW_DATA_DBL \
#      RRTMG_SW_DATA \
#      RRTMG_SW_DATA_DBL \
#      SOILPARM.TBL \
#      tr49t67 \
#      tr49t85 \
#      tr67t85 \
#      urban_param.tbl \
#      URBPARM.TBL \
#      VEGPARM.TBL \
#      '
#   foreach f ( $wrf_aux_file_list )
#      if ( -e $f ) then
#         ln -sf ${CENTRALDIR}/$f .
#         echo looking for  ${CENTRALDIR}/$f
#      else
#         echo $f not found.
#      endif
#   end
   ln -sf ${CENTRALDIR}/CAM_ABS_DATA .
   ln -sf ${CENTRALDIR}/CAM_AEROPT_DATA .
   ln -sf ${CENTRALDIR}/ETAMPNEW_DATA .
   ln -sf ${CENTRALDIR}/ETAMPNEW_DATA_DBL .
   ln -sf ${CENTRALDIR}/GENPARM.TBL .
   ln -sf ${CENTRALDIR}/grib2map.tbl .
   ln -sf ${CENTRALDIR}/gribmap.txt .
   ln -sf ${CENTRALDIR}/LANDUSE.TBL .
   ln -sf ${CENTRALDIR}/ozone.formatted .
   ln -sf ${CENTRALDIR}/ozone_lat.formatted .
   ln -sf ${CENTRALDIR}/ozone_plev.formatted .
   ln -sf ${CENTRALDIR}/RRTM_DATA .
   ln -sf ${CENTRALDIR}/RRTM_DATA_DBL .
   ln -sf ${CENTRALDIR}/RRTMG_LW_DATA .
   ln -sf ${CENTRALDIR}/RRTMG_LW_DATA_DBL .
   ln -sf ${CENTRALDIR}/RRTMG_SW_DATA .
   ln -sf ${CENTRALDIR}/RRTMG_SW_DATA_DBL .
   ln -sf ${CENTRALDIR}/SOILPARM.TBL .
   ln -sf ${CENTRALDIR}/tr49t67 .
   ln -sf ${CENTRALDIR}/tr49t85 .
   ln -sf ${CENTRALDIR}/tr67t85 .
   ln -sf ${CENTRALDIR}/URBPARM.TBL .
   ln -sf ${CENTRALDIR}/VEGPARM.TBL .

   # nfile is required when using mpi to run wrf.exe
   # nfile is machine specific.
   
   hostname >! nfile
   hostname >>! nfile
   ###ln -s  ${CENTRALDIR}/nfile$element nfile

   # Provides auxilliary info not available from DART state vector
   ${COPY} ${CENTRALDIR}/wrfinput_d0? .
   
   # if a mean state ic file exists convert it to a wrfinput_mean netcdf file
   if (  -e ${CENTRALDIR}/assim_model_state_ic_mean ) then
      ln -sf ${CENTRALDIR}/assim_model_state_ic_mean dart_wrf_vector
      ${CENTRALDIR}/dart_to_wrf >&! out.dart_to_wrf_mean
      ${COPY} wrfinput_d01 wrfinput_mean
      ${REMOVE} wrf.info dart_wrf_vector
   endif
   
   # ICs for this wrf run; Convert DART file to wrfinput netcdf file
   ${MOVE} ${CENTRALDIR}/$input_file dart_wrf_vector 
   
   ${CENTRALDIR}/dart_to_wrf >&! out.dart_to_wrf
   
   ${REMOVE} dart_wrf_vector
   
   # if append_lsm_data exists, run it.
   if ( -e ${CENTRALDIR}/append_lsm_data ) then
      ln -sf ${CENTRALDIR}/LSM/lsm_data_${element}.nc lsm_data.nc
      ${CENTRALDIR}/append_lsm_data
   endif

   # The program dart_tf_wrf has created the file wrf.info.
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
   # Set infl = 0. to avoid further deterioration of the ensemble member.
   
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
      if ( $BLOWN > 0 ) then
         set infl = 0.0
      endif
   endif
   
   # Find all BC's file available and sort them with "keys".
   # NOTE: this needs a fix for the idealized wrf case in which there are no
   # boundary files (also same for global wrf).  right now some of the
   # commands below give errors, which are ok to ignore in the idealized case
   # but it is not good form to generate spurious error messages.
   
   #--1st, check if LBCs are "specified" (in which case wrfbdy files are req'd)
   set SPEC_BC = `grep specified ${CENTRALDIR}/namelist.wrf.tmp | grep true | cat | wc -l`
   
   if ($SPEC_BC > 0) then
#      echo 'Specified BC in ADVANCE MODEL STATE'
      if ( -e ${CENTRALDIR}/da_wrfvar.exe ) then
#        echo 'Found da_wrfvar.exe'
        ls ${CENTRALDIR}/WRF/wrfbdy_*_mean* >! bdy.list
      else
      ls ${CENTRALDIR}/WRF/wrfbdy_*_$element >! bdy.list
      endif
      echo ${CENTRALDIR}/WRF/wrfbdy_d01_mean_ >! str.name
      sed 's/\//\\\//g' < str.name >! str.name2
      set STRNAME = `cat str.name2`
      set COMMAND = s/`echo ${STRNAME}`//
      sed $COMMAND < bdy.list >! bdy.list2
      sed 's/_/ /g' < bdy.list2 >! bdy.list
      set num_files = `cat bdy.list | wc -l`
      set items = `cat bdy.list`
      
      set ifile = 1
      set iday = 1
      set isec = 2
      if ( -e keys ) ${REMOVE} keys
      while ( $ifile <= $num_files )
         set key = `echo "$items[$iday] * 86400 + $items[$isec]" | bc`
         echo $key >>! keys
         set ifile = `expr $ifile \+ 1`
         set iday = `expr $iday \+ 2`
         set isec = `expr $isec \+ 2`
      end
      set keys = `sort keys`
  
   else
      echo ${CENTRALDIR}/WRF/wrfbdy_${targdays}_${targsecs}_$element >! bdy.list
   
   echo ${CENTRALDIR}/WRF/wrfbdy_ >! str.name
   sed 's/\//\\\//g' < str.name >! str.name2
   set STRNAME = `cat str.name2`
   set COMMAND = s/`echo ${STRNAME}`//
   
   sed $COMMAND < bdy.list >! bdy.list2
   sed 's/_/ /g' < bdy.list2 >! bdy.list
   set num_files = `cat bdy.list | wc -l`
   set items = `cat bdy.list`
   set ifile = 1
   set iday = 1
   set isec = 2
   if ( -e keys ) ${REMOVE} keys
   while ( $ifile <= $num_files )
      set key = `echo "$items[$iday] * 86400 + $items[$isec]" | bc`
      echo $key >>! keys
      set ifile = `expr $ifile \+ 1`
      set iday = `expr $iday \+ 3`
      set isec = `expr $isec \+ 3`
   end
   set keys = `sort keys`

   endif
   set keys = `sort keys`
  

 
   set cal_date    = `head -3 wrf.info | tail -1`
   set START_YEAR  = $cal_date[1]
   set START_MONTH = $cal_date[2]
   set START_DAY   = $cal_date[3]
   set START_HOUR  = $cal_date[4]
   set START_MIN   = $cal_date[5]
   set START_SEC   = $cal_date[6]
   
   set END_YEAR    = $cal_date[1]
   set END_MONTH   = $cal_date[2]
   set END_DAY     = $cal_date[3]
   set END_HOUR    = $cal_date[4]
   set END_MIN     = $cal_date[5]
   set END_SEC     = $cal_date[6]
   
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
      if ($ifile < $num_files ) then
         set ifile = `expr $ifile \+ 1`
      else
         echo No boundary file available to move beyond
         echo ${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:${START_MIN}:${START_SEC}
         exit
      endif
   end
   
   ###############################################################
   # Advance the model with new BC until target time is reached. #
   ###############################################################
   
   while ( $wrfkey < $targkey )
   
      set iday = `echo "$keys[$ifile] / 86400" | bc`
      set isec = `echo "$keys[$ifile] % 86400" | bc`
   
      # Copy the boundary condition file to the temp directory.
   
      if ( -e ${CENTRALDIR}/da_wrfvar.exe ) then
         ${COPY} ${CENTRALDIR}/WRF/wrfbdy_d01_mean_${iday}_${isec} wrfbdy_d01
      else
      ${COPY} ${CENTRALDIR}/WRF/wrfbdy_${iday}_${isec}_$element wrfbdy_d01
      endif
   
      if ( $targkey > $keys[$ifile] ) then
         set INTERVAL_SS = `echo "$keys[$ifile] - $wrfkey" | bc`
      else
         set INTERVAL_SS = `echo "$targkey - $wrfkey" | bc`
      endif
      set RUN_HOURS   = `expr $INTERVAL_SS \/ 3600`
      set REMAIN      = `expr $INTERVAL_SS \% 3600`
      set RUN_MINUTES = `expr $REMAIN \/ 60`
      set RUN_SECONDS = `expr $REMAIN \% 60`
      set INTERVAL_MIN = `expr $INTERVAL_SS \/ 60`
   
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
      set END_STRING = ${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:${END_MIN}:${END_SEC}
   
   # nreps should be the same as the number of domains you have, but
   # needs to be <= max_dom.  set max_dom up if you get an error.
   set nreps = ${MY_NUM_DOMAINS}

   #-----------------------------------------------------------------------
   # Update time control entries in the WRF namelist.input:
   #-----------------------------------------------------------------------
   
   ${COPY} ${CENTRALDIR}/namelist.wrf.tmp  namelist.input.hold
   ${REMOVE} script1.sed
   cat >! script1.sed << EOF
    /run_hours/c\
    run_hours                  = ${RUN_HOURS}
    /run_minutes/c\
    run_minutes                = ${RUN_MINUTES}
    /run_seconds/c\
    run_seconds                = ${RUN_SECONDS}
    /start_year/c\
    start_year                 = ${nreps}*${START_YEAR},
    /start_month/c\
    start_month                = ${nreps}*${START_MONTH},
    /start_day/c\
    start_day                  = ${nreps}*${START_DAY},
    /start_hour/c\
    start_hour                 = ${nreps}*${START_HOUR},
    /start_minute/c\
    start_minute               = ${nreps}*${START_MIN},
    /start_second/c\
    start_second               = ${nreps}*${START_SEC},
    /end_year/c\
    end_year                   = ${nreps}*${END_YEAR},
    /end_month/c\
    end_month                  = ${nreps}*${END_MONTH},
    /end_day/c\
    end_day                    = ${nreps}*${END_DAY},
    /end_hour/c\
    end_hour                   = ${nreps}*${END_HOUR},
    /end_minute/c\
    end_minute                 = ${nreps}*${END_MIN},
    /end_second/c\
    end_second                 = ${nreps}*${END_SEC},
# set history interval equal to run interval to make sure you have
# a wrfoutput file at the end of the run interval
    /history_interval/c\
    history_interval           = ${nreps}*${INTERVAL_MIN},
#  dart_to_wrf is expecting only a single time per file
    /frames_per_outfile/c\
    frames_per_outfile         = ${nreps}*1,
    /max_dom/c\
    max_dom                    = ${nreps},
EOF
   
       ${MOVE} namelist.input.hold namelist.input.tmp
       sed -f script1.sed namelist.input.tmp >! namelist.input.hold

       ${COPY} ${CENTRALDIR}/namelist.3dvar.input namelist.input
   
      # Update boundary conditions.
      # WARNING: da_wrfvar.exe will only work correctly if running WRF V3.1 or later!
      # If it is found in the central dir, use it to regnerate perturbed boundary files
      # Otherwise, do the original call to update_wrf_bc
      if ( -e ${CENTRALDIR}/da_wrfvar.exe ) then
         
#         echo 'da_wrfvar.exe found, using it to regenerate boundary conditions'
         ln -sf ${CENTRALDIR}/be.dat .
         set pscale = `head -1 ${CENTRALDIR}/bc_pert_scale | tail -1`
         set autoc  = `head -2 ${CENTRALDIR}/bc_pert_scale | tail -1`
         @ iseed2 = $element * 10000         
         set this_date     = ${START_YEAR}${START_MONTH}${START_DAY}${START_HOUR}
         echo this date is ${this_date}
         set da_window_min = `${CENTRALDIR}/da_advance_time.exe $this_date -3 -w`
         set da_window_max = `${CENTRALDIR}/da_advance_time.exe $this_date  3 -w`
         echo da_windows ${da_window_min} ${da_window_max}
         ${REMOVE} script.sed                
         cat >! script.sed << EOF
   /analysis_date/c\
   analysis_date = \'${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:${START_MIN}:${START_SEC}.0000\',
   /as1/c\
   as1 = ${pscale}, 2.0, 1.5,
   /as2/c\
   as2 = ${pscale}, 2.0, 1.5,
   /as3/c\
   as3 = ${pscale}, 2.0, 1.5,
   /as4/c\
   as4 = ${pscale}, 2.0, 1.5,
   /as5/c\
   as5 = ${pscale}, 2.0, 1.5,
   /var_scaling1/c\
   var_scaling1 = ${pscale},
   /var_scaling2/c\
   var_scaling2 = ${pscale},
   /var_scaling3/c\
   var_scaling3 = ${pscale},
   /var_scaling4/c\
   var_scaling4 = ${pscale},
   /var_scaling5/c\
   var_scaling5 = ${pscale},
   /seed_array1/c\
   seed_array1 = ${END_YEAR}${END_MONTH}${END_DAY}${END_HOUR},
   /seed_array2/c\
   seed_array2 = $iseed2,
   /time_window_min/c\
    time_window_min = "${da_window_min}.0000",
   /time_window_max/c\
    time_window_max = "${da_window_max}.0000",
   /start_year/c\
   start_year = ${START_YEAR},
   /start_month/c\
   start_month = ${START_MONTH},
   /start_day/c\
   start_day = ${START_DAY},
   /start_hour/c\
   start_hour = ${START_HOUR},
   /end_year/c\
   end_year = ${END_YEAR},
   /end_month/c\
   end_month = ${END_MONTH},
   /end_day/c\
   end_day = ${END_DAY},
   /end_hour/c\
   end_hour = ${END_HOUR},
EOF
      
         ${MOVE} namelist.input namelist.input.tmp
         sed -f script.sed namelist.input.tmp >! namelist.input
#         echo here is the first guess file name: ${CENTRALDIR}/WRF/wrfinput_d01_mean_${targdays}_${targsecs}
         ${COPY} ${CENTRALDIR}/WRF/wrfinput_d01_mean_${targdays}_${targsecs} ./fg
         ${CENTRALDIR}/da_wrfvar.exe >>&! out.wrfvar
         if ( -e rsl.out.0000 ) cat rsl.out.0000 >> out.wrfvar

         ${MOVE} wrfvar_output wrfinput_next
         ln -sf wrfinput_d01 wrfinput_this
         ln -sf wrfbdy_d01 wrfbdy_this

         # if wrfinput_mean file found, rename it
         if ( -e wrfinput_mean ) then
            ${MOVE} wrfinput_mean   wrfinput_this_mean
            ${MOVE} fg              wrfinput_next_mean
         endif

         echo $autoc | ${CENTRALDIR}/pert_wrf_bc >&! out.pert_wrf_bc
         ${REMOVE} wrfinput_this wrfinput_next wrfbdy_this
         if ( -e wrfinput_this_mean ) ${REMOVE} wrfinput_this_mean wrfinput_next_mean

      else  # Update boundary conditions

         echo 'Using update_wrf_bc to update boundary conditions'
         echo $infl | ${CENTRALDIR}/update_wrf_bc >&! out.update_wrf_bc

      endif
   
      # if you are running with nests, here is one way to get the info
      # into the wrf namelist.
      if ( -e ${CENTRALDIR}/tc_domain_info ) then

         set ndomains    = `head -1 ${CENTRALDIR}/tc_domain_info | tail -1`
         set i_start_str = `head -2 ${CENTRALDIR}/tc_domain_info | tail -1`
         set j_start_str = `head -3 ${CENTRALDIR}/tc_domain_info | tail -1`
         set nx_string   = `head -4 ${CENTRALDIR}/tc_domain_info | tail -1`
         set ny_string   = `head -5 ${CENTRALDIR}/tc_domain_info | tail -1`

         ${REMOVE} script.sed
         cat >! script.sed << EOF
         /max_dom/c\
         max_dom                             = ${ndomains},
         /e_we/c\
         e_we                                = ${nx_string},
         /e_sn/c\
         e_sn                                = ${ny_string},
         /i_parent_start/c\
         i_parent_start                      = ${i_start_str},
         /j_parent_start/c\
         j_parent_start                      = ${j_start_str},
EOF

         ${MOVE} namelist.input namelist.input.tmp
         sed -f script.sed namelist.input.tmp >! namelist.input

      endif

      if ( -e rsl.out.integration ) then
         ${REMOVE} rsl.*
      endif
   
      # run WRF here
      ${MOVE} namelist.input namelist.input.3dvar
      ${MOVE} namelist.input.hold namelist.input
      ${ADV_MOD_COMMAND} >>&! rsl.out.integration
      if ( -e rsl.out.0000 ) cat rsl.out.0000 >> rsl.out.integration
      ${COPY} rsl.out.integration $WRFOUTDIR/wrf.out_${targdays}_${targsecs}_${element}  
#      sleep 1
   
      set SUCCESS = `grep "wrf: SUCCESS COMPLETE WRF" rsl.* | cat | wc -l`
      if ($SUCCESS == 0) then
         echo $element >>! ${CENTRALDIR}/blown_${targdays}_${targsecs}.out
         echo "Model failure! Check file " ${CENTRALDIR}/blown_${targdays}_${targsecs}.out
         echo "for a list of failed elements, and check here for the individual output files:"
         echo " ${CENTRALDIR}/wrf.out_${targdays}_${targsecs}_${element}  "
         exit -1
      endif

      if ( -e ${CENTRALDIR}/extract ) then
         if ( $element == 1 ) then
            ls wrfout_d0${MY_NUM_DOMAINS}_* >! wrfout.list
            if ( -e ${CENTRALDIR}/psfc.nc ) then
               ${COPY} ${CENTRALDIR}/psfc.nc .
            endif
            echo `cat wrfout.list | wc -l` | ${CENTRALDIR}/extract
            ${MOVE} psfc.nc ${CENTRALDIR}/.
         endif
      endif
   
      set dn = 1
      while ( $dn <= $MY_NUM_DOMAINS )
         if ( $element <= $save_elements ) then
         ${COPY} wrfout_d0${dn}_${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:${END_MIN}:${END_SEC} $WRFOUTDIR/wrfout_d0${dn}_${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:${END_MIN}:${END_SEC}_$element
         endif
         ${MOVE} wrfout_d0${dn}_${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:${END_MIN}:${END_SEC} wrfinput_d0${dn}
         set dn = `expr $dn \+ 1`
      end
   
      ${REMOVE} wrfout*
   
      set START_YEAR  = $END_YEAR
      set START_MONTH = $END_MONTH
      set START_DAY   = $END_DAY
      set START_HOUR  = $END_HOUR
      set START_MIN   = $END_MIN
      set START_SEC   = $END_SEC
      set wrfkey = $keys[$ifile]
      set ifile = `expr $ifile \+ 1`
   
   end
   
   ##############################################
   # At this point, the target time is reached. #
   ##############################################
   
   # create new input to DART (taken from "wrfinput")
   ${CENTRALDIR}/wrf_to_dart >&! out.wrf_to_dart
   
   ${MOVE} dart_wrf_vector ${CENTRALDIR}/$output_file

   # and now repeat the entire process for any other ensemble member that
   # needs to be advanced by this task.
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
\rm -rf $control_file


