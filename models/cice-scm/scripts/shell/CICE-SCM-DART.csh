#!/bin/csh

set datea     = ${1}
#set datea     = 2011011300
set CASENAME='test11'
set NENS=80 
set NTRUTH=64 
set start_time=`date +%s`
set sim_length = 24
set storage_loc = /Users/bitz/models/ICEPACK_DART_RUNS
#set dart_loc = /Users/bitz/Dropbox/models/DART/models/cice-scm/work
set dart_loc = /Users/bitz/Dropbox/models/DARTmw/models/cice-scm2/work

set obs_dir = /Users/bitz/Dropbox/models/ICEPACK_DART_RUNS/OBS_SEQ_FILES/ITD11


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


###
echo "Starting Post Assimilation Processes for $datea"
###
set year = `echo $datea | cut -b1-4`
set month = `echo $datea | cut -b5-6`
set day = `echo $datea | cut -b7-8`
set hour = `echo $datea | cut -b9-10`
###
set fore_date = `echo $datea +24 | /Users/bitz/Dropbox/models/DART/models/cice/work/advance_time`
set fyear = `echo $fore_date | cut -b1-4`
set fmonth = `echo $fore_date | cut -b5-6`
set fday = `echo $fore_date | cut -b7-8`
set fhour = `echo $fore_date | cut -b9-10`
###
set past_date = `echo $datea -24 | /Users/bitz/Dropbox/models/DART/models/cice/work/advance_time`
###
set RUN_ASSIM = True
set RUN_POST_ASSIM = True

if ( ${RUN_ASSIM} == True) then
  echo "Running Assimilation..."
endif
if ( ${RUN_POST_ASSIM} == True) then
  echo "Running Post Assimilation processes, including icepack, ..."
endif

set run_dir = `pwd`
cd ${run_dir}

# do some prep work to improve speed
set mem = 1
while ($mem <= 80)
  set inst_string = `printf %04d $mem`
  cd mem$inst_string
  if (! -d analyses) then
    mkdir analyses
  endif
  cd ../
  @ mem++
end

if ( ${RUN_ASSIM} == True) then
  set n = 1
  while ($n <= ${NENS})
    set icnum = `echo $n + 10000 | bc | cut -b2-5`
    if (-e mem$icnum/restart/iced.${year}-${month}-${day}-00000.nc) then
      cd mem$icnum
      if ($n == ${NTRUTH}) then
        #echo "mv iced.${year}-${month}-${day}-00000.nc truth_data.nc"
        mv restart/iced.${year}-${month}-${day}-00000.nc truth_data.nc
      else
        #echo "mv iced.${year}-${month}-${day}-00000.nc dart_restart.nc"
        mv restart/iced.${year}-${month}-${day}-00000.nc dart_restart.nc
      endif
      cd ../
    endif
    @ n++
  end
  echo "Adding vars to dart_restart.nc"
  python ./addvars_to_restart.py
  rm -f obs_seq.out
  ln -s ${obs_dir}/obs_seq.out_${year}${month}${day} obs_seq.out
  # PRIOR INFLATION FILES LINKED IN
  #ln -s /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${past_date}/output_priorinf_mean.nc input_priorinf_mean.nc
  #ln -s /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${past_date}/output_priorinf_sd.nc input_priorinf_sd.nc
   
  echo "Running the filter..."
  ${dart_loc}/filter >& filter.out
  # on mac the stdio is buffered so the filter.out file is never complete
  # however dart_log.out is appended with most stuff that should go to filter.out and it is flushed
  # grep 'Finished' filter.out
  tail dart_log.out | grep 'Finished' 
  if ( $status != 0) then
    echo "Filter did not finish running"
    touch ${run_dir}/STOP
    exit
  else
    rm input_priorinf_*.nc
    echo "Filter is done running, modify restart_state.nc."
    python ./reconstitute_restart.py
  endif
endif
##
if ( ${RUN_POST_ASSIM} == True) then
# this is an example of syntax    if !(-d mem${inst_string}) mkdir -p mem${inst_string}
#  mkdir /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}
#  mkdir /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/analyses
#  mkdir /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts
#  mv filter.out /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/
#  mv regression_data.txt obs_inc_data.txt /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/
  ############
#  cp input_*.nc output_*.nc obs_seq.final /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/
####### PART I archive the dart statistics
  if (! -d dart_io) then
    mkdir dart_io
  endif
  #if (! -e dart_io/input_mean.${year}-${month}-${day}-00000.nc) then
    ncks -O -d ni,2,2 input_mean.nc dart_io/input_mean.${year}-${month}-${day}-00000.nc
    ncks -O -d ni,2,2 input_sd.nc dart_io/input_sd.${year}-${month}-${day}-00000.nc
    ncks -O -d ni,2,2 output_mean.nc dart_io/output_mean.${year}-${month}-${day}-00000.nc
    ncks -O -d ni,2,2 output_sd.nc dart_io/output_sd.${year}-${month}-${day}-00000.nc
    cp -f obs_seq.final dart_io/obs_seq.final_${year}${month}${day} 
  #endif

####### PART II run dart_to_cice and archive restarts
echo "Running dart_to_cice on each ensemble"
  set mem = 1
  while ($mem <= 80)
    set inst_string = `printf %04d $mem`
#    set icnum = `echo $mem + 10000 | bc | cut -b2-5`
    cd mem$inst_string
    if ($mem == 64) then
      #echo "mv truth_data.nc restart/iced.${year}-${month}-${day}-00000.nc"
      mv truth_data.nc restart/iced.${year}-${month}-${day}-00000.nc
    else
      ncks -O -d ni,2,2 dart_restart.nc ./analyses/iced.orig.${year}-${month}-${day}-00000.nc
      ncks -O -d ni,2,2 restart_state.nc ./analyses/restart_state.${year}-${month}-${day}-00000.nc
      cp ../input.nml .
      rm -f update_state.out
      #echo "/Users/bitz/Dropbox/models/DART/models/cice-scm/work/dart_to_cice >& update_state.out"
      ${dart_loc}/dart_to_cice >& update_state.out
      # on mac the stdio is buffered so the filter.out file is never complete
      # however dart_log.out is appended with most stuff that should go to filter.out and it is flushed
      # grep 'Finished' update_state.out
#      tail dart_log.out | grep 'Finished' 
#      if ( $status != 0) then
#        echo "Filter did not finish running"
#        touch ${run_dir}/STOP
#        exit
#      endif
      ncks -O -d ni,2,2 dart_restart.nc ./analyses/iced.updated.${year}-${month}-${day}-00000.nc
      #echo "mv dart_restart.nc restart/iced.${year}-${month}-${day}-00000.nc"
      mv dart_restart.nc restart/iced.${year}-${month}-${day}-00000.nc
    endif
    cd ../
    @ mem++
  end

####### PART III run icepack
echo "Running icepack on each ensemble"
  set mem = 1
  while ($mem <= 80)
    set inst_string = `printf %04d $mem`
#    set icnum = `echo $mem + 10000 | bc | cut -b2-5`
    cd mem$inst_string

   #rm icepack_in
   set restart_flag = .true.
   set startup_flag = .false.
   cat >! script.sed << EOF
    /RESTART_FILENAME/c\
    ice_ic         = 'restart/iced.${year}-${month}-${day}-00000.nc'
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
#      echo "./icepack >& icepack.out"
      ./icepack >& icepack.out
#    grep 'ICEPACK COMPLETED SUCCESSFULLY' icepack.out
#    if ($status != 0) then
#      echo "ICEPACK DID NOT FINISH....STOP!"
#      touch ${run_dir}/STOP
#      exit
#    endif   
#    if (! -e restart/iced.${fyear}-${fmonth}-${fday}-00000.nc) then
#      echo "ICEPACK DID NOT FINISH....STOP!"
#      touch ${run_dir}/STOP
#      exit
#    endif
    ########################
    cd ../
    @ mem++
  end
  rm -f mem*/restart_state.nc
  echo "Done with post assimilation step!"
endif
echo "Done cycling for $datea"
set end_time=`date +%s`
echo "execution time was `expr $end_time - $start_time` s."




