#!/bin/csh

# run a free run with initial conditions from the spinup restart
# which is iced.2012-01-01-00000.nc even though I cycled 20 times with year 2011 forcing
# the first year will be restart = .true. and always with ice_ic = iced.2012-01-01-00000.nc
# however because year_init is 2011, the date in icepack should be set to 2011-01-01
# and the first year will rewrite over iced.2012-01-01-00000.nc

# must have runtype_startup = 'T' for the first year


set storage_loc = /Users/bitz/models/ICEPACK_DART_RUNS
set sim_length = 24

set rsno_params  = ( -1.1662 -1.5390 -1.5601 -1.3936 -1.2739 -1.7979 -1.6640 -1.9532 -1.7008 -1.4012 \
-1.2438 -1.6094 -1.3161 -1.5490 -1.4130 -1.7756 -1.6207 -1.3761 -1.5467 -1.8103 \
-1.1176 -1.1707 -1.8719 -1.1091 -1.5732 -1.9218 -1.7375 -1.1296 -1.1492 -1.5762 \
-1.3516 -1.3200 -1.0940 -1.2564 -1.8513 -1.8130 -1.6986 -1.1607 -1.5495 -1.6511 \
-1.8117 -1.4098 -1.8445 -1.5602 -1.2460 -1.4585 -1.4374 -1.7378 -1.7524 -1.6781 \
-1.3135 -1.6909 -1.8998 -1.8430 -1.3014 -1.9527 -1.5665 -1.5061 -1.9653 -1.3262 \
-1.5849 -1.8442 -1.5210 -1.7941 -1.9308 -1.1359 -1.5627 -1.5851 -1.0821 -1.2352 \
-1.1719 -1.7196 -1.5995 -1.9995 -1.4933 -1.3038 -1.1022 -1.4972 -1.5062 -1.2231 )

#########################
set ksno_params = ( 0.2418 0.2797 0.2871 0.2892 0.2274 0.3015 0.2640 0.2647 0.2345 0.3129 \
0.2914 0.2982 0.2993 0.3193 0.2781 0.3000 0.2967 0.2550 0.2988 0.2175 \
0.3119 0.3119 0.2560 0.3004 0.2719 0.2886 0.2267 0.2838 0.2139 0.3038 \
0.2696 0.2630 0.2690 0.2325 0.3110 0.2864 0.3021 0.2178 0.3172 0.2401 \
0.2910 0.2851 0.2354 0.2511 0.2847 0.2890 0.2520 0.2577 0.3165 0.2512 \
0.2104 0.3177 0.2373 0.2767 0.2604 0.2353 0.2367 0.2286 0.2736 0.2422 \
0.3087 0.2890 0.3072 0.2542 0.2809 0.2893 0.2250 0.2195 0.2424 0.3190 \
0.2977 0.2629 0.2587 0.3003 0.2807 0.2950 0.2111 0.3045 0.2447 0.2760 )

set mem = 1
while ($mem <= 80)
  set inst_string = `printf %04d $mem`
  echo "running member $mem free run 2011-2020"
  mkdir mem${inst_string}
  mkdir mem${inst_string}/history
  mkdir mem${inst_string}/restart
  mkdir mem${inst_string}/save
  ln -s /Users/bitz/Dropbox/models/ICEPACK_DART_RUNS/TEST_SETUP/icepack mem${inst_string}/
  cp /Users/bitz/Dropbox/models/ICEPACK_DART_RUNS/spinup/mem${inst_string}/restart/iced.2012-01-01-00000.nc mem${inst_string}/restart
  cd mem${inst_string}
  set restart_flag = .true.
  set startup_flag = .true.    
  cat >! script.sed << EOF
    /RESTART_FILENAME/c\
    ice_ic         = 'restart/iced.2012-01-01-00000.nc'
    /SIM_LEN/c\
    npt              = $sim_length
    /RESTART_LOGIC/c\
    restart           = $restart_flag
    /STARTUP_LOGIC/c\
    runtype_startup   = $startup_flag
    /KSNOW/c\
    ksno              = ${ksno_params[$mem]}
    /RSNOW/c\
    R_snw             = ${rsno_params[$mem]}
    /MEMA/c\
    atm_data_file   = 'ATM_FORCING_${inst_string}.txt'
    /MEMO/c\
    ocn_data_file   = 'OCN_FORCING_${inst_string}.txt'  
EOF
  sed -f script.sed ../icepack_in.template >! icepack_in

  rm -f icepack.out script.sed
  ./icepack >& icepack.out
  grep 'ICEPACK COMPLETED SUCCESSFULLY' icepack.out
  if ($status != 0) then
    echo "ICEPACK DID NOT FINISH....STOP!"
    exit
  endif   
  if (! -e restart/iced.2011-01-02-00000.nc) then
    echo "ICEPACK DID NOT FINISH....STOP!"
    exit
  endif
#  ncks -d time,0,23,24 -d ni,2,2 history/icepack.h.20110101.nc save/icepack.h.20110101-00.nc
#  ncks -d time,23,23,24 -d ni,2,2 history/icepack.h.20110101.nc save/icepack.h.20110101-23.nc
#  rm history/icepack.h.20110101.nc 
  cd ../    
  @ mem = $mem + 1
end


