#!/bin/csh -f
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#----------------------------------------------------------------------
# This was written by Yongfei Zhang to explore running multiple
# instances of CICE with development version of CESM ... cesm1_5_beta06c
#----------------------------------------------------------------------

set myname = "yfzhang"
set casedir = "$WORK/cesmcases/"
set rundir  = "$SCRATCH/cesmruns/"

set num_instances = 30

setenv MYCASENAME  "Dtest_ens"$num_instances
setenv RES         T62_g16
setenv MACH        yellowstone 
setenv COMPSET     DTEST


cd $casedir

if ( -d $casedir/$MYCASENAME) then
    echo "case existed"
    echo "to delete the old case, use \
          rm -rf $casedir/$MYCASENAME" ;exit
endif

set scriptdir = "/glade/scratch/bitz/darttest/cesm1_5_beta06c/cime/scripts/"

${scriptdir}/create_newcase  -res  $RES \
                             -mach $MACH \
                             -compset $COMPSET \
                             -case $MYCASENAME \
                             -project P93300065

cd $MYCASENAME

set stream_year_align = 2000
set stream_year_first = 2000
set stream_year_last  = 2003

./xmlchange -file env_build.xml   -id CESMSCRATCHROOT   -val "$rundir/"
./xmlchange -file env_run.xml     -id DIN_LOC_ROOT      -val "$SCRATCH/inputdata_cam"
./xmlchange -file env_run.xml     -id RUNDIR            -val "$rundir/$MYCASENAME/run"
./xmlchange -file env_run.xml     -id STOP_OPTION       -val ndays
./xmlchange -file env_run.xml     -id STOP_N            -val 1
./xmlchange -file env_run.xml     -id RESUBMIT          -val 1

./xmlchange -file env_run.xml     -id DATM_MODE         -val CPLHIST3HrWx
./xmlchange -file env_run.xml     -id DATM_CPLHIST_YR_START -val $stream_year_first
./xmlchange -file env_run.xml     -id DATM_CPLHIST_YR_END   -val $stream_year_last
./xmlchange -file env_run.xml     -id DATM_CPLHIST_YR_ALIGN -val $stream_year_align
./xmlchange -file env_run.xml     -id RUN_STARTDATE     -val ${stream_year_first}-01-01

./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val $num_instances
./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val $num_instances
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val $num_instances
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 60
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 60
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 60
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 60
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 60
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 60
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 60
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 60


#./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1
#./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1
#./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1
#./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1
#./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1
#./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val 1
#./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val 1
#./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val 1

#./xmlchange -file env_mach_pes.xml -id TOTALPES            -val 64
#./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE  -val 16

./case.setup
echo "case setup finished"

@ inst = 1

while ( $inst <= $num_instances )

    set inst_string = `printf %04d $inst`

    set inst_2d = `printf %02d $inst`

    set fname = "user_nl_datm"_${inst_string}

    echo "streams  = 'datm.streams.txt.CAM4.multiyear.Solar_$inst_string             $stream_year_align $stream_year_first $stream_year_last'," >> ${fname}
    echo "           'datm.streams.txt.CAM4.multiyear.Precip_$inst_string            $stream_year_align $stream_year_first $stream_year_last'," >> ${fname}
    echo "           'datm.streams.txt.CAM4.multiyear.Other_$inst_string              $stream_year_align $stream_year_first $stream_year_last'"  >> ${fname}
    echo "           'datm.streams.txt.presaero.clim_2000_$inst_string 1 1 1'"  >> ${fname}
    echo "vectors  = 'u:v' "     >> ${fname}
    echo "mapmask  = 'nomask', " >> ${fname}
    echo "           'nomask', " >> ${fname}
    echo "           'nomask', " >> ${fname}
    echo "           'nomask'  " >> ${fname}
    echo "tintalgo = 'coszen', " >> ${fname}
    echo "           'nearest'," >> ${fname}
    echo "           'linear', " >> ${fname}
    echo "           'linear'  " >> ${fname}


    # Create stream files for each ensemble member
    cp $WORK/Stream_Tempolate/datm.streams.txt.CAM4.multiyear.Precip \
       datm.streams.txt.CAM4.multiyear.Precip_${inst_string}
    cp $WORK/Stream_Tempolate/datm.streams.txt.CAM4.multiyear.Solar \
       datm.streams.txt.CAM4.multiyear.Solar_${inst_string}
    cp $WORK/Stream_Tempolate/datm.streams.txt.CAM4.multiyear.Other \
       datm.streams.txt.CAM4.multiyear.Other_${inst_string}

    foreach FNAME (datm.streams.txt.CAM4.multiyear.*)
      sed s/CAM_DATM.cpl/CAM_DATM-${inst_2d}.cpl/ $FNAME >! temp
      mv temp $FNAME
    end

    # Creat user_nl_cice for each ensemble member
    cp datm.streams.txt.CAM4.multiyear.* ${rundir}/$MYCASENAME/run
 
    set fname = "user_nl_cice_${inst_string}"
#    echo "ice_ic    = '/glade/scratch/yfzhang/inputdata_cam/ice/cice/sp2000_ens1.cice.r.2043-01-01-00000.nc' " >>$fname
    echo "histfreq_n =  1,1,1,1,1  "                             >> $fname
    echo "histfreq   = 'd','m','x','x','x' "                               >> $fname
    echo "f_sst = 'dmxxx' " >> $fname
    echo "f_sss = 'dmxxx' " >> $fname
    echo "f_frzmlt = 'dmxxx' " >> $fname
    echo "f_frz_onset = 'dmxxx' " >>$fname
    echo "f_aicen = 'dmxxx' " >>$fname

    @ inst = $inst + 1
end

./case.build

exit 0


