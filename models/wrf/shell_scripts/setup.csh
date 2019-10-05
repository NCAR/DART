#!/bin/csh
#
# This script sets up the proper directory structure for running the WRF/DART tutorial.
# It also places the required DART, WRF, and WRFDA executables and files in the proper location.
#
# >> ./setup.csh param.csh
#

set paramfile = `readlink -f ${1}` # Get absolute path for param.csh from command line arg
source $paramfile
set COPY = 'cp -rL'

# Create directories
foreach dir ( ${RUN_DIR} ${TEMPLATE_DIR} ${OBSPROC_DIR} ${OUTPUT_DIR} ${ICBC_DIR} \
              ${POST_STAGE_DIR} ${OBS_DIAG_DIR} ${PERTS_DIR} ${SHELL_SCRIPTS_DIR} \
              ${RUN_DIR}/WRFIN ${RUN_DIR}/WRFOUT ${RUN_DIR}/WRF_RUN )
    if ( ! -d $dir ) mkdir -p $dir
end

# Put DART executables in proper place
foreach exe ( advance_time filter pert_wrf_bc obs_diag obs_sequence_tool \
              obs_seq_to_netcdf fill_inflation_restart wrf_dart_obs_preprocess )
    if ( -e ${DART_DIR}/models/wrf/work/$exe ) then
       ${COPY} ${DART_DIR}/models/wrf/work/$exe ${RUN_DIR} || exit 1
    else
       echo "ERROR: ${DART_DIR}/models/wrf/work/$exe not found; cannot be copied to ${RUN_DIR}"
       echo "Check DART build."
       exit
    endif
end

if ( -e ${DART_DIR}/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc ) then
   ${COPY} ${DART_DIR}/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc ${RUN_DIR} || exit 2
    else
       echo "ERROR: ${DART_DIR}/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc does not exit"
       echo "Check DART directory."
       exit
    endif
endif

# Put WRF executables and supporting files in proper place
${COPY} ${WRF_DM_SRC_DIR}/run/* ${RUN_DIR}/WRF_RUN  || exit 3
if ( ! -e ${RUN_DIR}/WRF_RUN/wrf.exe ) then
    echo "ERROR: wrf.exe not copied into ${RUN_DIR}/WRF_RUN"
    echo "Check WRF build"
endif
${REMOVE} ${RUN_DIR}/WRF_RUN/namelist.input  || exit 4

if ( -e ${WRF_SL_SRC_DIR}/main/real.exe ) then
   ${COPY} ${WRF_SL_SRC_DIR}/main/real.exe ${RUN_DIR}/WRF_RUN/real.serial.exe || exit 5
else
   echo "ERROR: ${WRF_SL_SRC_DIR}/main/real.exe not found; cannot be copied to ${RUN_DIR}/WRF_RUN/"
   echo "Check WRF build."
   exit
endif

# Put WRFDA executables and supporting files in proper place
if ( -e ${VAR_SRC_DIR}/var/build/da_wrfvar.exe ) then
  ${COPY} ${VAR_SRC_DIR}/var/build/da_wrfvar.exe ${RUN_DIR}/WRF_RUN  || exit 6
else
   echo "ERROR: ${VAR_SRC_DIR}/var/build/da_wrfvar.exe not found; cannot be copied to ${RUN_DIR}/WRF_RUN/"
   echo "Check WRFDA build."
   exit
endif 
if ( -e ${VAR_SRC_DIR}/var/run/be.dat.cv3 ) then
  ${COPY} ${VAR_SRC_DIR}/var/run/be.dat.cv3 ${RUN_DIR}/WRF_RUN/be.dat || exit 7
else
   echo "ERROR: ${VAR_SRC_DIR}/var/run/be.dat.cv3 not found; cannot be copied to ${RUN_DIR}/WRF_RUN/"
   echo "Check WRFDA build."
   exit
endif 

# Put scripts in proper place 
${COPY} ${SHELL_SCRIPTS_DIR}/add_bank_perts.ncl ${RUN_DIR} || exit 8
${COPY} ${SHELL_SCRIPTS_DIR}/new_advance_model.csh ${RUN_DIR} || exit 9

exit


