#!/bin/csh -f
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#-----------------------------------------------------------------------
# Script init_ens_remote4.csh
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
setenv START_DATE_ASSIM 2002110200            # Start time of period.
setenv NCYCLE 40                              # Number of assimilation cycles.
setenv FCST_RANGE 6                           # Forecast range (hours).
setenv INTERVAL 6                             # Interval between analyses (hours)
setenv WRF_DT 600                             # Model timestep (seconds)

setenv AVN_DIR     /ocotillo1/wrfdev/AVN      # Global analysis directory
setenv SRC_DIR     /ocotillo1/${user}         # Location of codes:
setenv DAT_DIR     ${SRC_DIR}/GEN_INIT_ENS    # Scratch data space.

setenv WRFSI_DIR_SRC   ${SRC_DIR}/wrfsi_20020326  # WRF SI.
setenv WRFSI_DIR   ${WRFSI_DIR_SRC}           # WRF SI.

#setenv WRFSI_DIR /usr/local/wrfsi/wrfsi_20020328 # WRFSI.
setenv WRFSI_DAT_DIR /usr/local/wrfsi/SI_GEOG    # WRFSI input data.
setenv WRF_DIR     ${SRC_DIR}/WRFV1.3            # WRF

# End of user modifications.

setenv INSTALLROOT $WRFSI_DIR_SRC
setenv FCST_RANGE_SEC `expr $FCST_RANGE \* 3600`
setenv NUM_TIMESTEPS `expr $FCST_RANGE_SEC \/ ${WRF_DT}`
setenv OUT_FREQ $NUM_TIMESTEPS

set startnode = 2
set endnode = 10

set ES = 81

echo "${ES}" > ens_size

rm ${WRF_DIR}/test/em_real/run_real_*.out

set seconds = $FCST_RANGE_SEC
set days = 0

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
	    msrcp mss:/DSS/DS083.2/data/$AVN_FILE $AVN_DIR/.
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
   set WRFSI_DIR = /var/tmp/wrfsi_20020326_${NC}
   (rsh -n node$inode "rm -rf ${WRFSI_DIR} >& /dev/null" )
   rcp -r $WRFSI_DIR_SRC node${inode}:${WRFSI_DIR}
   set INSTALLROOT = $WRFSI_DIR

   set workdir = `pwd`
   (rsh -n node$inode "cd $workdir; run_wrfsi_remote.csh $START_DATE $FCST_RANGE $INTERVAL $DAT_DIR_MEM $WRFSI_DAT_DIR $INSTALLROOT $WRFSI_DIR >& preprocess_${NC}" ) &

    set START_DATE = `advance_cymdh ${START_DATE} 24`   # Go to next member (day)

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

   (rsh -n node$inode "rm -rf /var/tmp/GEN_INIT_ENS_${NC}" )

#--------------------------------------------------------------------
# 3) Second stage of preprocessing into WRF format (runs WRF real.exe). 
#--------------------------------------------------------------------

   rm ${WRF_DIR}/test/em_real/real_input_em*
   rcp -r node${inode}:/var/tmp/wrfsi_20020326_${NC}/data/siprd ${WRF_DIR}/test/em_real/.
   mv ${WRF_DIR}/test/em_real/siprd/real_input_em* ${WRF_DIR}/test/em_real/.
   rm -r ${WRF_DIR}/test/em_real/siprd

   run_wrfreal_remote.csh $NC $START_DATE $FCST_RANGE $INTERVAL $WRF_DT $NUM_TIMESTEPS $OUT_FREQ $WRF_DIR

   (rsh -n node$inode "rm -rf /var/tmp/wrfsi_20020326_${NC}")
   mv ${WRF_DIR}/test/em_real/wrfinput_d01 ${DAT_DIR}/wrfinput_${NC}
   mv ${WRF_DIR}/test/em_real/wrfbdy_d01 ${DAT_DIR}/wrfbdy_${NC}
   mv ${WRF_DIR}/test/em_real/namelist.input ${DAT_DIR}/namelist.input_${days}_${seconds}_${NC}

   set START_DATE = `advance_cymdh ${START_DATE} 24`   # Go to next member (day)

   @ NC ++

   if ($inode == $endnode) then
      set inode = $startnode
   else
      @ inode ++
   endif

end

if ( $ICYC == 1 ) then
   cp ${DAT_DIR}/wrfinput_21 ${DAT_DIR}/wrfinput_mean
endif
cp ${DAT_DIR}/wrfbdy_21 ${DAT_DIR}/wrfbdy_mean

ensemble_init < ens_size > out.ensemble_init_${days}_${seconds}

mv ${DAT_DIR}/wrfbdy_mean ${DAT_DIR}/wrfbdy_mean_${days}_${seconds}

# Prepare to do the ensemble initial conditions file.
if ( $ICYC == 1 ) then
   echo ".false." >  input_wrf_to_dart
   echo "0" > time.dat
   echo "0" >> time.dat

   rm -f wrf_ics
endif

set NC = 1
# Loop over the ensemble members
while ( $NC <= $ES )
   if ( $ICYC == 1 ) then
#---------------------------------------------------
# Convert wrfinput (netcdf) files into dart readable
#---------------------------------------------------
      mv ${DAT_DIR}/wrfinput_${NC} wrfinput
      
# create new input to DART (taken from "wrfinput")
      dart_tf_wrf < input_wrf_to_dart > out.wrf_to_dart

      cat dart_wrf_vector >> wrf_ics

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

set seconds = `expr $seconds \+ $FCST_RANGE_SEC`
if ( $seconds >= 86400) then
   set seconds = `expr $seconds \- 86400`
   set days = `expr $days \+ 1`
endif

end

exit (0)
