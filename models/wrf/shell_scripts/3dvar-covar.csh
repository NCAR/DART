#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Purpose: Given a first guess ensemble mean, generate ensemble members
#          from wrf/3dvar covariances.

set echo

#--------------------------------------------
# 0) Set up various environment variables:
#--------------------------------------------

set ES = 84

set seconds = 0
set days = 146827

setenv MY_NUM_DOMAINS 1

# End of user modifications.

echo $seconds $days > wrf.info

set SEED1 = 1
set dn = 1
while ( $dn <= $MY_NUM_DOMAINS )

   cp wrfinput_d0${dn}_mean_${days}_${seconds} wrf_3dvar_input

   set NC = 1
# Loop over the ensemble members
   while ( $NC <= $ES )

@ SEED2 = ${SEED1} * 100

rm -f script.sed
cat > script.sed << EOF
 s/SEED1/${SEED1}/
 s/SEED2/${SEED2}/
EOF

 sed -f script.sed \
    namelist.3dvar.template > namelist.3dvar

      ./da_3dvar.exe >& da_3dvar.out_${dn}_${NC}
#      ./wrfvar.exe >& wrfvar.out_${dn}_${NC}

      mv wrf_3dvar_output wrfinput_d0${dn}_${NC}

      @ NC ++

      @ SEED1 ++

   end

   @ dn ++

end

rm -f filter_ics

set NC = 1
# Loop over the ensemble members
while ( $NC <= $ES )

   set dn = 1
   while ( $dn <= $MY_NUM_DOMAINS )
      mv wrfinput_d0${dn}_${NC} wrfinput_d0${dn}
      @ dn ++
   end

#---------------------------------------------------
# Convert wrfinput (netcdf) files into dart readable
#---------------------------------------------------

# create new input to DART (taken from "wrfinput_d0x")
   wrf_to_dart >& out.wrf_to_dart

   cat dart_wrf_vector >> filter_ics

   rm dart_wrf_vector

   @ NC ++

end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

