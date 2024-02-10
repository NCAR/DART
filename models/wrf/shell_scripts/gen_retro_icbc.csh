#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#==================================================================
#BSUB -J gen_retro_icbc
#BSUB -o gen_retro_icbc.%J.log
#BSUB -P 25000077 # if using LSF, change this
#BSUB -W 0:38
#BSUB -q regular
#BSUB -n 64
#BSUB -x
#BSUB -R "span[ptile=64]"

#PBS -N gen_retro_icbc
#PBS -A 25000077 # if using PBS, change this
#PBS -l walltime=00:38
#PBS -q regular               
#PBS -o gen_retro_icbc.out                    
#PBS -j oe                              
#PBS -k eod                              
#PBS -l select=5:ncpus=60:mpiprocs=60
#PBS -V                        

echo "gen_retro_icbc.csh is running in `pwd`"

########################################################################
#
#   gen_retro_icbc.csh - shell script that generates the
#                                     necessary wrfinput_d01 and
#                                     wrfbdy_d01 files for running
#                                     a real-time analysis system.
#
#     created May 2009, Ryan Torn, U. Albany
#
# This creates   output/${date}/wrfbdy_d01_{days}_{seconds}_mean
#                output/${date}/wrfinput_d01_{days}_{time_step1}_mean
#                output/${date}/wrfinput_d01_{days}_{time_step2}_mean
########################################################################

set datea     = 2017042700
set datefnl   = 2017042712 # set this appropriately #%%%#
set paramfile = /glade/derecho/scratch/USERNAME/WORK_DIR/scripts/param.csh   # set this appropriately #%%%#

source $paramfile

# The geo_*.nc files should already be in the ${ICBC_DIR}/*/ directories.
# ${LINK} ${GEO_FILES_DIR}/geo_*.nc .

mkdir -p ${ICBC_DIR}/metgrid
${LINK} ${WPS_SRC_DIR}/metgrid/METGRID.TBL ${ICBC_DIR}/metgrid/METGRID.TBL

while ( 1 == 1 )
   echo "Entering gen_retro_icbc.csh for $datea"

   if ( ! -d ${OUTPUT_DIR}/${datea} )  mkdir -p ${OUTPUT_DIR}/${datea}

   cd ${ICBC_DIR}
   ${LINK} ${RUN_DIR}/input.nml input.nml
   ${REMOVE} gfs*pgrb2* *grib2

   #  prepare to run WPS ungrib and metgrid
   set start_date = `echo $datea 0 -w | ${DART_DIR}/models/wrf/work/advance_time`
   set end_date   = `echo $datea 6 -w | ${DART_DIR}/models/wrf/work/advance_time`
   echo $start_date

   ${REMOVE} script.sed
   ${REMOVE} namelist.wps
   cat >! script.sed << EOF
 /start_date/c\
 start_date = 2*'${start_date}',
 /end_date/c\
 end_date   = 2*'${end_date}',
EOF

   # build grib file names - may need to change for other data sources.

   if ( $GRIB_SRC != 'GFS' ) then 
      echo "gen_retro_icbc.csh: GRIB_SRC is set to $GRIB_SRC"
      echo "gen_retro_icbc.csh: There are some assumptions about using 'GFS'."
      echo "If you want to use something else, you will need to change the"
      echo "values of gribfile_a, gribfile_b"
      exit 2
   endif

   set gribfile_a = ${GRIB_DATA_DIR}/${datea}/gfs_ds084.1/gfs.0p25.${datea}.f000.grib2
   set gribfile_b = ${GRIB_DATA_DIR}/${datea}/gfs_ds084.1/gfs.0p25.${datea}.f006.grib2
   ${LINK} $gribfile_a GRIBFILE.AAA
   ${LINK} $gribfile_b GRIBFILE.AAB

   sed -f script.sed ${TEMPLATE_DIR}/namelist.wps.template >! namelist.wps
   ${LINK} ${WPS_SRC_DIR}/ungrib/Variable_Tables/Vtable.${GRIB_SRC} Vtable

   ${REMOVE}                    output.ungrib.exe.${GRIB_SRC}
   ${WPS_SRC_DIR}/ungrib.exe >& output.ungrib.exe.${GRIB_SRC}

   ${REMOVE}                     output.metgrid.exe
   ${WPS_SRC_DIR}/metgrid.exe >& output.metgrid.exe

   ${LINK} ${WPS_SRC_DIR}/met_em.d01.* .
   set datef  =  `echo $datea $ASSIM_INT_HOURS | ${DART_DIR}/models/wrf/work/advance_time`
   set gdatef = (`echo $datef 0 -g             | ${DART_DIR}/models/wrf/work/advance_time`)
   set hh     =  `echo $datea | cut -b9-10`

   #  Run real.exe twice, once to get first time wrfinput_d0? and wrfbdy_d01,
   #  then again to get second time wrfinput_d0? file
   set n = 1
   while ( $n <= 2 )

      echo "RUNNING REAL, STEP $n"
      echo " "

      if ( $n == 1 ) then
         set date1      = $datea
         set date2      = $datef
         set fcst_hours = $ASSIM_INT_HOURS
      else
         set date1      = $datef
         set date2      = $datef
         set fcst_hours = 0
      endif

      set yyyy1 = `echo $date1 | cut -c 1-4`
      set mm1   = `echo $date1 | cut -c 5-6`
      set dd1   = `echo $date1 | cut -c 7-8`
      set hh1   = `echo $date1 | cut -c 9-10`
      set yyyy2 = `echo $date2 | cut -c 1-4`
      set mm2   = `echo $date2 | cut -c 5-6`
      set dd2   = `echo $date2 | cut -c 7-8`
      set hh2   = `echo $date2 | cut -c 9-10`

      ${REMOVE} namelist.input script.sed
      cat >! script.sed << EOF
    /run_hours/c\
    run_hours                  = ${fcst_hours},
    /run_minutes/c\
    run_minutes                = 0,
    /run_seconds/c\
    run_seconds                = 0,
    /start_year/c\
    start_year                 = ${yyyy1}, ${yyyy1},
    /start_month/c\
    start_month                = ${mm1}, ${mm1},
    /start_day/c\
    start_day                  = ${dd1}, ${dd1},
    /start_hour/c\
    start_hour                 = ${hh1}, ${hh1},
    /start_minute/c\
    start_minute               = 00, 00,
    /start_second/c\
    start_second               = 00, 00,
    /end_year/c\
    end_year                   = ${yyyy2}, ${yyyy2},
    /end_month/c\
    end_month                  = ${mm2}, ${mm2},
    /end_day/c\
    end_day                    = ${dd2}, ${dd2},
    /end_hour/c\
    end_hour                   = ${hh2}, ${hh2},
    /end_minute/c\
    end_minute                 = 00, 00,
    /end_second/c\
    end_second                 = 00, 00,
EOF

      sed -f script.sed ${TEMPLATE_DIR}/namelist.input.meso >! namelist.input
     #${REMOVE} out.real.exe
     #${RUN_DIR}/WRF_RUN/real.serial.exe >& out.real.exe
     #if ( -e rsl.out.0000 )  cat rsl.out.0000 >> out.real.exe

      rm script.sed real_done rsl.*
      echo "2i\"                                       >! script.sed
      echo "#======================================\"  >> script.sed
      echo "#PBS -N run_real\"                         >> script.sed
      echo "#PBS -A ${COMPUTER_CHARGE_ACCOUNT}\"       >> script.sed
      echo "#PBS -l walltime=00:05:00\"                >> script.sed
      echo "#PBS -q ${ADVANCE_QUEUE}\"                 >> script.sed
      echo "#PBS -l job_priority=${ADVANCE_PRIORITY}\" >> script.sed
      echo "#PBS -o run_real.out\"                     >> script.sed
      echo "#PBS -j oe\"                               >> script.sed
      echo "#PBS -k eod\"                              >> script.sed
      echo "#PBS -l select=1:ncpus=128:mpiprocs=128\"  >> script.sed
      echo "#PBS -V\"                                  >> script.sed
      echo "#======================================\"  >> script.sed
      echo "\"                                         >> script.sed
      echo ""                                          >> script.sed
      echo 's%${1}%'"${paramfile}%g"                   >> script.sed
      sed -f script.sed ${SHELL_SCRIPTS_DIR}/real.csh  >! real.csh

      qsub real.csh

      # need to look for sometihng to know when this job is done
      while ( ! -e ${ICBC_DIR}/real_done )
          sleep 15
      end
      cat rsl.out.0000 >> out.real.exe

      #  move output files to storage
      set gdate = (`echo $date1 0 -g | ${DART_DIR}/models/wrf/work/advance_time`)
      ${MOVE} wrfinput_d01 ${OUTPUT_DIR}/${datea}/wrfinput_d01_${gdate[1]}_${gdate[2]}_mean
      if ( $n == 1 ) ${MOVE} wrfbdy_d01 ${OUTPUT_DIR}/${datea}/wrfbdy_d01_${gdatef[1]}_${gdatef[2]}_mean

      @ n++

   end

   # move to next time, or exit if final time is reached
   if ( $datea == $datefnl) then
      echo "Reached the final date "
      echo "Script exiting normally"
      exit 0
   endif
   set datea  = `echo $datea $ASSIM_INT_HOURS | ${DART_DIR}/models/wrf/work/advance_time`
   echo "starting next time: $datea"
end

exit 0

