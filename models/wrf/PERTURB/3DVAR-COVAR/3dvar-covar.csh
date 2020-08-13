#!/bin/csh -f
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Purpose: Given a first guess ensemble mean, generate ensemble members
#          from wrf/3dvar covariances.
#
# What you need in the local directory:
#
# ensemble mean:                wrfinput_d01_mean, wrfinput_d02_mean, ...
# wrf/3dvar executable:         da_3dvar.exe
# background error statistics:  be.cv_2, be.cv_3
# wrf/3dvar namelists:          namelist.input, namelist.3dvar.template
# dart_tf_wrf converter:        dart_tf_wrf
# dart_tf_wrf namelist:         input.nml
#
# At the beginning of this script, specify
#
# 1. the ensemble size (ES);
# 2. Valid time of the ensemble (seconds, days);
# 3. Number of wrf domains (MY_NUM_DOMAINS).
#
# Also, namelist.input and input.nml should be edited according to your needs.
#
# The namelist.3dvar.template is used by the present script. DO NOT change the lines
#
# seed_array1    = SEED1,
# seed_array2    = SEED2/
#
# in the record11 section.
#-----------------------------------------------------------------------

set echo

#--------------------------------------------
# Set up various environment variables:
#--------------------------------------------

set ES = 40

set seconds = 0
set days = 146827

setenv MY_NUM_DOMAINS 1

# End of user modifications.

rm -f fort.33
ln -s be.cv_3 fort.33

echo $seconds $days > wrf.info

set SEED1 = 1
set dn = 1
while ( $dn <= $MY_NUM_DOMAINS )

   cp wrfinput_d0${dn}_mean wrf_3dvar_input

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
   echo ".false." | dart_tf_wrf >& out.wrf_to_dart

   cat dart_wrf_vector >> filter_ics

   rm dart_wrf_vector

   @ NC ++

end

exit 0


