#!/bin/csh -f
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#-----------------------------------------------------------------------
# Script init_ens_ocotillo.csh
#
# Purpose: Given start date etc, runs whole WRF system from global 
#          analysis (NCEP AVN/FNL):
#
#          1) Prepare global analysis (NCEP/FNL - 1deg res.).
#          2) Preprocess into WRF format (runs SI and real).
#
#-----------------------------------------------------------------------

#--------------------------------------------
# 0) Set up various environment variables:
#--------------------------------------------

setenv RPC_UNSUPPORTED_NETIFS eth0

setenv DOMAIN CONUS90                         # Domain name.
setenv RUN_ID TEST                            # Character string.
#setenv START_DATE_ASSIM 2002110200            # Start time of period.
setenv START_DATE_ASSIM 2003010100            # Start time of period.
setenv NCYCLE 4                              # Number of assimilation cycles.
setenv FCST_RANGE 6                           # Forecast range (hours).
setenv INTERVAL 6                             # Interval between analyses (hours)
setenv WRF_DT 600                             # Model timestep (seconds)

setenv AVN_DIR     /ocotillo1/wrfdev/AVN      # Global analysis directory

setenv S_AVAIL_DATE 1999122500
setenv E_AVAIL_DATE 2004060700
setenv BACK_YEAR       5000000

setenv SRC_DIR     /ocotillo1/${user}         # Location of codes:
setenv DAT_DIR     ${SRC_DIR}/GEN_ENS         # Scratch data space.

setenv WRFSI_DIR_SRC  ${SRC_DIR}/WRFV2/wrfsi  # WRF SI.
setenv WRFSI_DIR   ${WRFSI_DIR_SRC}           # WRF SI.

#setenv WRFSI_DIR /usr/local/wrfsi/wrfsi_20020328 # WRFSI.
setenv WRF_DIR     ${SRC_DIR}/WRFV2            # WRF

set startnode = 1
set endnode = 10

set ES = 80

 setenv WEST_EAST_GRIDS	  45
 setenv SOUTH_NORTH_GRIDS 45
 setenv VERTICAL_GRIDS	  28
 setenv GRID_DISTANCE	  200000

set seconds = 0
set days = 146827

echo $seconds $days > wrf.info

# End of user modifications.

setenv INSTALLROOT $WRFSI_DIR_SRC
setenv FCST_RANGE_SEC `expr $FCST_RANGE \* 3600`
setenv OUT_FREQ `expr $FCST_RANGE_SEC \/ ${WRF_DT}`

set seconds = `expr $seconds \+ $FCST_RANGE_SEC`
if ( $seconds >= 86400) then
   set seconds = `expr $seconds \- 86400`
   set days = `expr $days \+ 1`
endif

# Write the target time at the beginning of wrf.info

echo $seconds $days > temp
cat wrf.info >> temp
mv temp wrf.info

echo "${ES}" > ens_size

rm ${WRF_DIR}/test/em_real/run_real_*.out

set ICYC = 1
# Loop over cycles
while ( $ICYC <= $NCYCLE )

##########if ( $ICYC == 2 ) then

   set START_DATE = $START_DATE_ASSIM

   set inode = $startnode
   set cyclenode = 1
   set NC = 1
# Loop over the ensemble members
   while ( $NC <= $ES )

   set END_DATE = `advance_cymdh $START_DATE $FCST_RANGE`  # End time of forecast.
   set DAT_DIR_MEM = /var/tmp/GEN_INIT_ENS_${NC}

#---------------------------------------------------
# 1) Prepare global analysis (NCEP/FNL - 1deg res.).
#---------------------------------------------------

    echo "1) Prepare global analysis (NCEP/FNL - 1deg res.)."

    set DATE = $START_DATE
    (rsh -n node$inode "rm -rf ${DAT_DIR_MEM} >& /dev/null ; mkdir ${DAT_DIR_MEM} ; mkdir ${DAT_DIR_MEM}/AVN" )

    while ( $DATE <= $END_DATE )
	set YY = `echo $DATE | cut -c3-4`
	set MM = `echo $DATE | cut -c5-6`
	set DD = `echo $DATE | cut -c7-8`
	set HH = `echo $DATE | cut -c9-10`

	set AVN_FILE = fnl_${YY}${MM}${DD}_${HH}_00

	if ( -e ${AVN_DIR}/$AVN_FILE ) then
	    echo "   File $AVN_FILE exists in ${AVN_DIR}/"
	else 
	    echo "   Retrieving $AVN_FILE to $AVN_DIR"
#	    msrcp mss:/DSS/DS083.2/data/$AVN_FILE $AVN_DIR/.
            rsh -n bay "rm -f /mmmtmp/caya/migs/$AVN_FILE"
            rsh -n bay "msrcp mss:/DSS/DS083.2/data/$AVN_FILE /mmmtmp/caya/migs/. ; rcp /mmmtmp/caya/migs/$AVN_FILE ocotillo:$AVN_DIR/."
            rsh -n bay "rm -f /mmmtmp/caya/migs/$AVN_FILE"
	endif
        rcp ${AVN_DIR}/$AVN_FILE node${inode}:${DAT_DIR_MEM}/AVN/$AVN_FILE

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
   set WRFSI_DIR = /var/tmp/wrfsi_${NC}
   (rsh -n node$inode "rm -rf ${WRFSI_DIR} >& /dev/null" )
   rcp -r $WRFSI_DIR_SRC node${inode}:${WRFSI_DIR}
   set INSTALLROOT = $WRFSI_DIR

   set workdir = `pwd`
   (rsh -n node$inode "cd $workdir; ./run_wrfsi.csh $START_DATE $FCST_RANGE $INTERVAL $DAT_DIR_MEM $INSTALLROOT $WEST_EAST_GRIDS $SOUTH_NORTH_GRIDS $VERTICAL_GRIDS $GRID_DISTANCE >& preprocess_${NC}" ) &

    if (`expr ${START_DATE} \+ 1000000` <= $E_AVAIL_DATE) then
       set START_DATE = `expr ${START_DATE} \+ 1000000` # Go to next year
    else
       set START_DATE = `expr ${START_DATE} \- ${BACK_YEAR}`
       while (${START_DATE} < ${S_AVAIL_DATE})
          set START_DATE = `expr ${START_DATE} \+ 1000000` # Go to next year
       end
       set START_DATE = `advance_cymdh ${START_DATE} 72`     # Go to next member (+3days)
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

   rm ${WRF_DIR}/test/em_real/wrf_real_input_em*
   rcp -r node${inode}:/var/tmp/GEN_INIT_ENS_${NC}/data/siprd ${WRF_DIR}/test/em_real/.

########   (rsh -n node$inode "rm -rf /var/tmp/GEN_INIT_ENS_${NC}" )

   mv ${WRF_DIR}/test/em_real/siprd/wrf_real_input_em* ${WRF_DIR}/test/em_real/.
   rm -r ${WRF_DIR}/test/em_real/siprd

   ./run_wrfreal.csh $NC $START_DATE $FCST_RANGE $INTERVAL $WRF_DT $OUT_FREQ $WRF_DIR $WEST_EAST_GRIDS $SOUTH_NORTH_GRIDS $VERTICAL_GRIDS $GRID_DISTANCE

##############   (rsh -n node$inode "rm -rf /var/tmp/wrfsi_${NC}")

   mv ${WRF_DIR}/test/em_real/wrfinput_d01 ${DAT_DIR}/wrfinput_${NC}
   mv ${WRF_DIR}/test/em_real/wrfbdy_d01 ${DAT_DIR}/wrfbdy_${NC}
#   mv ${WRF_DIR}/test/em_real/namelist.input ${DAT_DIR}/namelist.input_${days}_${seconds}_${NC}
   rm ${WRF_DIR}/test/em_real/namelist.input

    if (`expr ${START_DATE} \+ 1000000` <= $E_AVAIL_DATE) then
       set START_DATE = `expr ${START_DATE} \+ 1000000` # Go to next year
    else
       set START_DATE = `expr ${START_DATE} \- ${BACK_YEAR}`
       while (${START_DATE} < ${S_AVAIL_DATE})
          set START_DATE = `expr ${START_DATE} \+ 1000000` # Go to next year
       end
       set START_DATE = `advance_cymdh ${START_DATE} 72`   # Go to next member (+3days)
    endif

   @ NC ++

   if ($inode == $endnode) then
      set inode = $startnode
   else
      @ inode ++
   endif

end

# Prepare to do the ensemble initial conditions file.
if ( $ICYC == 1 ) then
   cp ${DAT_DIR}/wrfinput_1 ${DAT_DIR}/wrfinput_mean
   rm -f filter_ics
endif
cp ${DAT_DIR}/wrfbdy_1 ${DAT_DIR}/wrfbdy_mean

ensemble_init < ens_size > out.ensemble_init_${days}_${seconds}

mv ${DAT_DIR}/wrfbdy_mean ${DAT_DIR}/wrfbdy_mean_${days}_${seconds}

set NC = 1
# Loop over the ensemble members
while ( $NC <= $ES )
   if ( $ICYC == 1 ) then
#---------------------------------------------------
# Convert wrfinput (netcdf) files into dart readable
#---------------------------------------------------
      mv ${DAT_DIR}/wrfinput_${NC} wrfinput

# create new input to DART (taken from "wrfinput")
      echo ".false." | dart_tf_wrf >& out.wrf_to_dart

      cat dart_wrf_vector >> filter_ics

      rm dart_wrf_vector
   else
      rm ${DAT_DIR}/wrfinput_${NC}
   endif
   mv ${DAT_DIR}/wrfbdy_${NC} ${DAT_DIR}/wrfbdy_${days}_${seconds}_${NC}
@ NC ++
end

##########endif

set START_DATE_ASSIM = `advance_cymdh $START_DATE_ASSIM $FCST_RANGE` # Advance to next cycle
@ ICYC ++

echo $seconds $days > wrf.info

set seconds = `expr $seconds \+ $FCST_RANGE_SEC`
if ( $seconds >= 86400) then
   set seconds = `expr $seconds \- 86400`
   set days = `expr $days \+ 1`
endif

# Write the target time at the beginning of wrf.info

echo $seconds $days > temp
cat wrf.info >> temp
mv temp wrf.info

end

exit (0)
