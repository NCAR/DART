#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Purpose: Given start date etc, runs whole WRF system from global 
#          analysis (NCEP AVN/FNL):
#
#          1) Prepare global analysis (NCEP/FNL - 1deg res.).
#          2) Preprocess into WRF format (runs SI and real).
#
#-----------------------------------------------------------------------

#      The GRID GENeration is the most expensive operation here.
#      Thus, it is advised to first run with
#
#      setenv RUN_GRIDGEN 1
#      NCYCLE = 1
#      ES = 1
#
#      This generate the GRID once and for all.
#      If the result is correct, then run with
#
#      setenv RUN_GRIDGEN 0
#
#      and increase NCYCLE and ES as desired.

#--------------------------------------------
# 0) Set up various environment variables:
#--------------------------------------------

setenv RUN_GRIDGEN 0
setenv RPC_UNSUPPORTED_NETIFS eth0

setenv CASE_NAME	     conus200
setenv START_DATE_ASSIM 2003010100            # Start time of period.

set ini_seconds = 0
set ini_days = 146827

setenv NCYCLE 4                             # Number of assimilation cycles.
setenv FCST_RANGE 6                           # Forecast range (hours).
setenv INTERVAL 6                             # Interval between analyses (hours)

setenv GRIB_DATA AVN                          # AVN = (NCEP/FNL - 1deg res.)
                                              # ETA = (GCIP NCEP Eta - ~40 km res.)

setenv DATASOURCE AVN

setenv DATA_DIR     /ocotillo6/caya/wrfdev/${GRIB_DATA}      # Global analysis directory

setenv S_AVAIL_DATE 1999122500
setenv E_AVAIL_DATE 2004060700
setenv BACK_YEAR       5000000

setenv DAT_DIR     `pwd`                      # Scratch data space.

#setenv WRF_DIR        /ocotillo1/${USER}/WRFV2.0.3.1/WRFV2   # WRF
setenv WRF_DIR        /ocotillo1/${USER}/WRFV2   # WRF
setenv WRFSI_DIR_SRC  ${WRF_DIR}/wrfsi        # WRF SI.

setenv VARTMP         /var/tmp                # Remote work directory

set startnode = 1
set endnode = 14

set ES = 1
set SCALE = 0.2

setenv WRF_DT 600                             # Model timestep (seconds)
setenv WEST_EAST_GRIDS	  45
setenv SOUTH_NORTH_GRIDS  45
setenv VERTICAL_GRIDS	  28
setenv GRID_DISTANCE	  200000
setenv MY_MOAD_KNOWN_LAT	"40.0"
setenv MY_MOAD_KNOWN_LON	"-98.0"
setenv MY_MOAD_STAND_LONS	"-98.0"
setenv MY_NUM_DOMAINS     1

setenv SF_SURFACE_PHYSICS 1

if ( $SF_SURFACE_PHYSICS == 1 ) then
   setenv NUM_SOIL_LAYERS 5
else
   setenv NUM_SOIL_LAYERS 4
endif

set LLI = (   1 39 39 )
set LLJ = (   1 36 36 )
set URI = (  45 78 78 )
set URJ = (  45 72 72 )

set grid_ratio = ( 1 3 3 )
set time_step_ratio = ( 1 3 3 )

set nextmem = 72                # Advance this many hours for the next member.

# End of user modifications.

set seconds = $ini_seconds
set days = $ini_days

set e_we = ( 0 0 0 )
set e_sn = ( 0 0 0 )
set dx = ( 0 0 0 )
set dt = ( 0 0 0 )

set dn = 1
while ( $dn <= $MY_NUM_DOMAINS )
   @ e_we[$dn] = ($URI[$dn] - $LLI[$dn]) * $grid_ratio[$dn] + 1
   @ e_sn[$dn] = ($URJ[$dn] - $LLJ[$dn]) * $grid_ratio[$dn] + 1
   if ( $dn == 1 ) then
      set dx[$dn] = $GRID_DISTANCE
      set dt[$dn] = $WRF_DT
      set dxn1 = $GRID_DISTANCE
      set dtn1 = $WRF_DT
   else
      @ dx[$dn] = $dxn1 / $grid_ratio[$dn]
      set dxn1 = $dx[$dn]
      @ dt[$dn] = $dtn1 / $time_step_ratio[$dn]
      set dtn1 = $dt[$dn]
   endif
   @ dn ++
end

echo $seconds $days > wrf.info

setenv INSTALLROOT $WRFSI_DIR_SRC
setenv FCST_RANGE_SEC `expr $FCST_RANGE \* 3600`
setenv OUT_FREQ `expr $FCST_RANGE_SEC \/ ${WRF_DT}`

set seconds = `expr $seconds \+ $FCST_RANGE_SEC`
if ( $seconds >= 86400) then
   set seconds = `expr $seconds \- 86400`
   set days = `expr $days \+ 1`
endif

# Write the target time at the end of wrf.info

echo $seconds $days >> wrf.info

echo "${ES}"     > ens.info
echo "${SCALE}" >> ens.info

rm ${WRF_DIR}/test/em_real/run_real_*.out

set ICYC = 1
# Loop over cycles
while ( $ICYC <= $NCYCLE )

   set START_DATE = $START_DATE_ASSIM

   set inode = $startnode
   set cyclenode = 1
   set NC = 1
# Loop over the ensemble members
   while ( $NC <= $ES )

   set END_DATE = `advance_cymdh $START_DATE $FCST_RANGE`  # End time of forecast.
   set DAT_DIR_MEM = ${VARTMP}/${user}_GEN_INIT_ENS_${NC}

#---------------------------------------------------
# 1) Prepare global analysis.
#---------------------------------------------------

   echo "1) Prepare global analysis from ${GRIB_DATA}."

   set DATE = $START_DATE
   (rsh -n node$inode "rm -rf ${DAT_DIR_MEM} >& /dev/null ; mkdir -p ${DAT_DIR_MEM}/${GRIB_DATA}" )
   while ( $DATE <= $END_DATE )

      set MM = `echo $DATE | cut -c5-6`
      set DD = `echo $DATE | cut -c7-8`
      set HH = `echo $DATE | cut -c9-10`

      if ($GRIB_DATA == AVN) then
	 set YY = `echo $DATE | cut -c3-4`
	 set GRIB_FILE = fnl_${YY}${MM}${DD}_${HH}_00
      else
	 set YY = `echo $DATE | cut -c1-4`
	 set GRIB_FILE = ${YY}${MM}${DD}${HH}.AWIP3D00.tm00
      endif

      if ( -e ${DATA_DIR}/$GRIB_FILE ) then
         echo "   File $GRIB_FILE exists in ${DATA_DIR}/"
      else 
         if ($GRIB_DATA == AVN) then
	    echo "   Retrieving $GRIB_FILE to $DATA_DIR"
#	     msrcp mss:/DSS/DS083.2/data/$GRIB_FILE $DATA_DIR/.
            rsh -n bay "rm -f /mmmtmp/${USER}/migs/$GRIB_FILE"
            rsh -n bay "msrcp mss:/DSS/DS083.2/data/$GRIB_FILE /mmmtmp/${USER}/migs/. ; rcp /mmmtmp/${USER}/migs/$GRIB_FILE ocotillo:$DATA_DIR/."
            rsh -n bay "rm -f /mmmtmp/${USER}/migs/$GRIB_FILE"
	 else 
	    echo "Please put $GRIB_FILE to $DATA_DIR first."
            exit
	 endif
      endif

      set DATE=`advance_cymdh ${DATE} ${INTERVAL}`

   end
   echo ""

#------------------------------------------------------------
# 2) First stage of preprocessing into WRF format (runs SI). 
#------------------------------------------------------------

   echo "2) First stage of preprocessing into WRF format (runs SI)."

# Clean out SI output files:
   rm ${WRFSI_DIR_SRC}/data/siprd/real_input_em* >& /dev/null 

   rm -f preprocess_${NC} >& /dev/null

   set workdir = `pwd`
   (rsh -n node$inode "cd $workdir; ./run_wrfsi.csh $START_DATE $FCST_RANGE $INTERVAL $DAT_DIR_MEM $INSTALLROOT $WEST_EAST_GRIDS $SOUTH_NORTH_GRIDS $VERTICAL_GRIDS $GRID_DISTANCE $MY_MOAD_KNOWN_LAT $MY_MOAD_KNOWN_LON $MY_MOAD_STAND_LONS $MY_NUM_DOMAINS $LLI[2] $LLJ[2] $URI[2] $URJ[2] $DATASOURCE $CASE_NAME $GRIB_DATA $RUN_GRIDGEN $DATA_DIR >& preprocess_${NC}" ) &

    if (`expr ${START_DATE} \+ 1000000` <= $E_AVAIL_DATE) then
       set START_DATE = `expr ${START_DATE} \+ 1000000`           # Go to next year
    else
       set START_DATE = `expr ${START_DATE} \- ${BACK_YEAR}`
       while (${START_DATE} < ${S_AVAIL_DATE})
          set START_DATE = `expr ${START_DATE} \+ 1000000`        # Go to next year
       end
       set START_DATE = `advance_cymdh ${START_DATE} ${nextmem}`  # Go to next member
    endif

   @ NC ++

   if ($inode == $endnode) then
      set inode = $startnode
      if ($cyclenode == 1) then
         @ cyclenode ++
      else
         set cyclenode = 1
         wait
      endif
   else
      @ inode ++
   endif

#   sleep 3

end

echo "Waiting for all processes to finish up..."
wait
sleep 10

set START_DATE = $START_DATE_ASSIM

set inode = $startnode
set NC = 1
# Loop over the ensemble members
while ( $NC <= $ES )

#--------------------------------------------------------------------
# 3) Second stage of preprocessing into WRF format (runs WRF real.exe). 
#--------------------------------------------------------------------

#-----------------------------------------------------------------------
# Create WRF namelist.input:
#-----------------------------------------------------------------------

source ./get_date_range.csh $START_DATE $FCST_RANGE

setenv INTERVAL_SS `expr $INTERVAL \* 3600`

rcp -r node${inode}:${VARTMP}/${user}_GEN_INIT_ENS_${NC}/data/siprd \
                    ${WRF_DIR}/test/em_real/.
########   (rsh -n node$inode "rm -rf ${VARTMP}/${user}_GEN_INIT_ENS_${NC}" )

set dn = 1
while ( $dn <= $MY_NUM_DOMAINS )

if ( $dx[$dn] < 5000 ) then
   setenv CU_PHYSICS 0
else
   setenv CU_PHYSICS 1
endif

cat >! ${WRF_DIR}/test/em_real/namelist.input << EOF
 &time_control
 run_days                   = 0,
 run_hours                  = ${FCST_RANGE},
 run_minutes                = 0,
 run_seconds                = 0,
 start_year                 = ${START_YEAR},
 start_month                = ${START_MONTH},
 start_day                  = ${START_DAY},
 start_hour                 = ${START_HOUR},
 start_minute               = 00,
 start_second               = 00,
 end_year                   = ${END_YEAR},
 end_month                  = ${END_MONTH},
 end_day                    = ${END_DAY},
 end_hour                   = ${END_HOUR},
 end_minute                 = 00,
 end_second                 = 00,
 interval_seconds           = ${INTERVAL_SS},
 input_from_file            = .true.,
 history_interval           = 360,
 frames_per_outfile         = 1000,
 restart                    = .false.,
 restart_interval           = 5000,
 io_form_history            = 2
 io_form_restart            = 2
 io_form_input              = 2
 io_form_boundary           = 2
 debug_level                = 0
 /

 &domains
 time_step                  = $dt[$dn],
 time_step_fract_num        = 0,
 time_step_fract_den        = 1,
 max_dom                    = 1,
 s_we                       = 1,
 e_we                       = $e_we[$dn],
 s_sn                       = 1,
 e_sn                       = $e_sn[$dn],
 s_vert                     = 1,
 e_vert                     = ${VERTICAL_GRIDS},
 dx                         = $dx[$dn],
 dy                         = $dx[$dn],
 grid_id                    = 1,
 parent_id                  = 0,
 i_parent_start             = 0,
 j_parent_start             = 0,
 parent_grid_ratio          = 1,
 parent_time_step_ratio     = 1,
 feedback                   = 1,
 smooth_option              = 1
 /

 &physics
 mp_physics                 = 3,
 ra_lw_physics              = 1,
 ra_sw_physics              = 1,
 radt                       = 30,
 sf_sfclay_physics          = 1,
 sf_surface_physics         = ${SF_SURFACE_PHYSICS},
 bl_pbl_physics             = 1,
 bldt                       = 0,
 cu_physics                 = ${CU_PHYSICS},
 cudt                       = 5,
 isfflx                     = 1,
 ifsnow                     = 0,
 icloud                     = 1,
 surface_input_source       = 1,
 num_soil_layers            = ${NUM_SOIL_LAYERS},
 maxiens                    = 1,
 maxens                     = 3,
 maxens2                    = 3,
 maxens3                    = 16,
 ensdim                     = 144,
 sst_update                 = 1
/

 &dynamics
 dyn_opt                    = 2,
 rk_ord                     = 3,
 w_damping                  = 0,
 diff_opt                   = 0,
 km_opt                     = 1,
 damp_opt                   = 0,
 zdamp                      = 5000.,
 dampcoef                   = 0.2,
 khdif                      = 0,
 kvdif                      = 0,
 smdiv                      = 0.1,
 emdiv                      = 0.01,
 epssm                      = 0.1,
 non_hydrostatic            = .true.,
 time_step_sound            = 4,
 h_mom_adv_order            = 5,
 v_mom_adv_order            = 3,
 h_sca_adv_order            = 5,
 v_sca_adv_order            = 3,
 /

 &bdy_control
 spec_bdy_width             = 5,
 spec_zone                  = 1,
 relax_zone                 = 4,
 specified                  = .true.,
 periodic_x                 = .false.,
 symmetric_xs               = .false.,
 symmetric_xe               = .false.,
 open_xs                    = .false.,
 open_xe                    = .false.,
 periodic_y                 = .false.,
 symmetric_ys               = .false.,
 symmetric_ye               = .false.,
 open_ys                    = .false.,
 open_ye                    = .false.,
 nested                     = .false.,
 /

 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
 /

EOF

   rm ${WRF_DIR}/test/em_real/wrf_real_input_em*

   mv -v ${WRF_DIR}/test/em_real/siprd/wrf_real_input_em.d0${dn}.* \
         ${WRF_DIR}/test/em_real/.
   if ($dn > 1) then
      mv -v ${WRF_DIR}/test/em_real/wrf_real_input_em.d0${dn}.${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00 \
            ${WRF_DIR}/test/em_real/wrf_real_input_em.d01.${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00
   endif

   ./run_wrfreal.csh $NC $WRF_DIR

   if ( $NC == 1 ) then
      mv -v ${WRF_DIR}/test/em_real/wrflowinp_d0${dn} \
            ${DAT_DIR}/wrflowinp_d0${dn}_${days}_${seconds}
   endif

   mv -v ${WRF_DIR}/test/em_real/wrfinput_d01 \
         ${DAT_DIR}/wrfinput_d0${dn}_${NC}
   if ($dn == 1) then
      mv -v ${WRF_DIR}/test/em_real/wrfbdy_d01 \
            ${DAT_DIR}/wrfbdy_${NC}
   endif

   @ dn ++

end

   rm -r ${WRF_DIR}/test/em_real/siprd
########   (rsh -n node$inode "rm -rf ${VARTMP}/${user}_wrfsi_${NC}")

   if (`expr ${START_DATE} \+ 1000000` <= $E_AVAIL_DATE) then
      set START_DATE = `expr ${START_DATE} \+ 1000000`            # Go to next year
   else
      set START_DATE = `expr ${START_DATE} \- ${BACK_YEAR}`
      while (${START_DATE} < ${S_AVAIL_DATE})
         set START_DATE = `expr ${START_DATE} \+ 1000000`         # Go to next year
      end
      set START_DATE = `advance_cymdh ${START_DATE} ${nextmem}`   # Go to next member
   endif

   @ NC ++

   if ($inode == $endnode) then
      set inode = $startnode
   else
      @ inode ++
   endif

end

# Save a copy of wrfinput for TSK, TMN, SST, VEGFRA, ALBBCK, IVGTYP
#set dn = 1
#while ( $dn <= $MY_NUM_DOMAINS )
#   cp -pv ${DAT_DIR}/wrfinput_d0${dn}_1 \
#          ${DAT_DIR}/wrfinput_d0${dn}.${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00
#   @ dn ++
#end

# Prepare to do the ensemble initial conditions file.
set dn = 1
while ( $dn <= $MY_NUM_DOMAINS )
   cp -pv ${DAT_DIR}/wrfinput_d0${dn}_1 \
          ${DAT_DIR}/wrfinput_d0${dn}_mean
   @ dn ++
end
if ( $ICYC == 1 ) then
   set dn = 1
   while ( $dn <= $MY_NUM_DOMAINS )
      cp -pv ${DAT_DIR}/wrfinput_d0${dn}_mean \
             ${DAT_DIR}/wrfinput_d0${dn}_mean_${ini_days}_${ini_seconds}
      @ dn ++
   end
   rm -f filter_ics
endif
cp -pv ${DAT_DIR}/wrfbdy_1 ${DAT_DIR}/wrfbdy_mean

ensemble_init < ens.info > out.ensemble_init_${days}_${seconds}

mv -v ${DAT_DIR}/wrfbdy_mean ${DAT_DIR}/wrfbdy_mean_${days}_${seconds}

set NC = 1
# Loop over the ensemble members
while ( $NC <= $ES )

   set dn = 1
   while ( $dn <= $MY_NUM_DOMAINS )
      mv -v ${DAT_DIR}/wrfinput_d0${dn}_${NC} wrfinput_d0${dn}
      @ dn ++
   end

   if ( $ICYC == 1 ) then
#---------------------------------------------------
# Convert wrfinput (netcdf) files into dart readable
#---------------------------------------------------

# create new input to DART (taken from "wrfinput_d0x")
      echo ".false." | dart_tf_wrf >& out.wrf_to_dart

      cat dart_wrf_vector >> filter_ics

      rm dart_wrf_vector

   endif

   mv -v ${DAT_DIR}/wrfbdy_${NC} \
         ${DAT_DIR}/wrfbdy_${days}_${seconds}_${NC}

   @ NC ++

end   # Loop over the ensemble members

set START_DATE_ASSIM = `advance_cymdh $START_DATE_ASSIM $FCST_RANGE` # Advance to next cycle

echo $seconds $days > wrf.info

set seconds = `expr $seconds \+ $FCST_RANGE_SEC`
if ( $seconds >= 86400) then
   set seconds = `expr $seconds \- 86400`
   set days = `expr $days \+ 1`
endif

# Write the target time at the end of wrf.info

echo $seconds $days >> wrf.info

@ ICYC ++

end   # Loop over cycles

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

