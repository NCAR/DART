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
       ${COPY} ${DART_DIR}/models/wrf/work/$exe ${RUN_DIR}
    else
       echo "ERROR: ${DART_DIR}/models/wrf/work/$exe not found; cannot be copied to ${RUN_DIR}"
       echo "Check DART build."
       exit
    endif
end

${COPY} ${DART_DIR}/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc ${RUN_DIR}

# Put WRF executables and supporting files in proper place
${COPY} ${WRF_DM_SRC_DIR}/run/* ${RUN_DIR}/WRF_RUN  ; ${REMOVE} ${RUN_DIR}/WRF_RUN/namelist.input #KRF remove namelist.input? it is overwritten later though
${COPY} ${WRF_SL_SRC_DIR}/main/real.exe ${RUN_DIR}/WRF_RUN/real.serial.exe

# Put WRFDA executables and supporting files in proper place
${COPY} ${VAR_SRC_DIR}/var/build/da_wrfvar.exe ${RUN_DIR}/WRF_RUN 
${COPY} ${VAR_SRC_DIR}/var/run/be.dat.cv3 ${RUN_DIR}/WRF_RUN/be.dat


# Put scripts in proper place #KRF Remove from init_ensemeble_var if we keep here
${COPY} ${SHELL_SCRIPTS_DIR}/add_bank_perts.ncl ${RUN_DIR}
${COPY} ${SHELL_SCRIPTS_DIR}/new_advance_model.csh ${RUN_DIR}


# KRF Right now the namelists and io lists are create/copied to rundir when executing
# init_ensemble_var.csh. Could do this in here in setup.csh, but they're case specific.

exit


