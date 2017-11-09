#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# The DART program 'filter' has the ability to take a single model state and 
# perturb it to create an ensemble of model states.
#
# The intent of the this script is to take a single coamps state and make N ensemble
# members using the DART program 'filter'. Each ensemble member will differ from the 
# parent in that noise will be added to each of the variables defined in the state.vars.
# The nature of the noise is defined in the filter namelist. It is not expected that 
# these files are valid for coamps, but they will enable testing of the assimilation.
#
# FIXME: Someone should take the logic from perturb_init.f90 and put it in 
# model_mod:pert_model_copies(). At present, pert_model_copies() uses a default routine.
#
# The steps are basically:
# 1) create an empty directory for the ensemble and copy all the bits needed to
#    create the ensemble into it.
# 2) run trans_coamps_to_dart to create a netCDF file - required by DART
# 3) create an observation sequence file (required by filter) by running pmo
# 4) create a bunch of output files that filter will update.
# 5) convert all those output files into coamps hdf5 restart files
# 6) clean up the garbage?
#
#===============================================================================
# This block of directives constitutes the preamble for the SLURM queuing system.
# the normal way to submit with slurm:  sbatch perturb_single.csh
# and to check on the job status     :  squeue --users=thoar
#
#SBATCH --job-name=ensemble
#SBATCH --output=ensemble%A.log
#SBATCH --error=ensemble%A.err
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH --error=ensemble%A.err
#SBATCH --output=ensemble%A.log
#===============================================================================

if ($?SLURM_JOB_ID) then
   set ORIGINALDIR = $SLURM_SUBMIT_DIR
   set     JOBNAME = $SLURM_JOB_NAME
   set       JOBID = $SLURM_JOBID
   set     MYQUEUE = $SLURM_JOB_PARTITION
   set      MYHOST = $SLURM_SUBMIT_HOST
   set    NODELIST = $SLURM_NODELIST
   set   LAUNCHCMD = "mpirun -np $SLURM_NTASKS -bind-to core"
else
   set ORIGINALDIR = `pwd`
   set     JOBNAME = ensemble
   set       JOBID = $$
   set     MYQUEUE = Interactive
   set      MYHOST = $host
   set   LAUNCHCMD = " "
endif

#==============================================================================
# Just an echo of job attributes
#==============================================================================

echo
echo "${JOBNAME} ($JOBID) submit directory ${ORIGINALDIR}"
echo "${JOBNAME} ($JOBID) submit      host ${MYHOST}"
echo "${JOBNAME} ($JOBID) running in queue ${MYQUEUE}"
echo "${JOBNAME} ($JOBID) running       on ${NODELIST}"
echo "${JOBNAME} ($JOBID) started at "`date`
echo

set ORIGIN = /home/nopp/COAMPS_hdf5_files
set DESTINATION = /home/thoar/COAMPS_hdf5_files/Ensemble2
set ENSEMBLE_SIZE = 4
set TEMPLATE = coamps_2013011000.hdf5

#==============================================================================
# Must convert a single state to a netcdf file for filter, perturb it,
# take the output netcdf files and update the hdf5 files ...
#==============================================================================

mkdir -p ${DESTINATION}
cd ${DESTINATION}
\rm -f dart_log.out dart_log.nml

ln -sf /home/thoar/DART/coamps/models/coamps_nest/work/trans_coamps_to_dart .   || exit 1
ln -sf /home/thoar/DART/coamps/models/coamps_nest/work/trans_dart_to_coamps .   || exit 1
ln -sf /home/thoar/DART/coamps/models/coamps_nest/work/perfect_model_obs    .   || exit 1
ln -sf /home/thoar/DART/coamps/models/coamps_nest/work/filter               .   || exit 1
cp     /home/thoar/DART/coamps/models/coamps_nest/work/input.nml            .   || exit 1
cp     /home/thoar/DART/coamps/models/coamps_nest/work/convert.nml          .   || exit 1
cp     /home/thoar/DART/coamps/models/coamps_nest/work/state.vars.full      state.vars   || exit 1
cp     /home/thoar/DART/coamps/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc .

# enforce the assumptions

ex input.nml <<ex_end
g;input_state_files ;s;= .*;= 'dart_vector.nc', 'dart_vector.nc';
g;input_state_file_list ;s;= .*;= '', '';
g;output_state_file_list ;s;= .*;= 'output_list_domain_1.txt', 'output_list_domain_1.txt';
g;ens_size ;s;= .*;= ${ENSEMBLE_SIZE};
g;num_output_obs_members ;s;= .*;= ${ENSEMBLE_SIZE};
g;num_output_state_members ;s;= .*;= 0;
g;output_mean ;s;= .*;= .FALSE.;
g;output_sd ;s;= .*;= .FALSE.;
g;perturb_from_single_instance ;s;= .*;= .TRUE.;
g;sampling_error_correction ;s;= .*;= .FALSE.;
wq
ex_end

# trans_coamps_to_dart creates a netCDF file (dart_vector.nc) from 
# an HDF5 file (coamps.hdf5).  DART reads the variable shapes, etc. from a 
# file called 'coamps.hdf5' but the actual input data files are specified 
# via namelist mechanisms.

cp   ${ORIGIN}/${TEMPLATE} .
ln -s          ${TEMPLATE} coamps.hdf5

./trans_coamps_to_dart  || exit 2

# create an observation sequence file (needed for filter)
# obs_seq.2obs.in has precisely 2 observations - one is identically
# part of the COAMPS state (at whatever index 2101 happens to be)
# and the other is a potential temperature on LEVEL 50.0 (wherever that is)
#
# The namelist for pmo reads the two nests, each from the same file: dart_vector.nc
# pmo is not set up to WRITE both nests to the same file when creating the output
# file from scratch, but will happily write to an existing file - so we copy
# the input file to an output file and let pmo overwrite it.

cp /home/thoar/DART/coamps/models/coamps_nest/work/obs_seq.2obs.in      obs_seq.in
cp dart_vector.nc perfect_output.nc

./perfect_model_obs || exit 3

# we need an output ensemble that will get updated with the perturbed states.

set instance = 1
while ( $instance <= ${ENSEMBLE_SIZE} )
   set OUTFILE = `printf dart_%03d_output.nc $instance`
   cp -v dart_vector.nc  ${OUTFILE}
   @ instance ++
end

#  create the list of output files for DART - one per domain. In this case
#  they are both the same since both domains exist in the same input file.
ls -1 dart_???_output.nc > output_list_domain_1.txt
ls -1 dart_???_output.nc > output_list_domain_2.txt

${LAUNCHCMD} ./filter || exit 4

# We now have an ensemble of netCDF files. We need an ensemble of HDF5 files.
# At present, the best way to do this is to copy the input file to a bunch
# of output files and let 'trans_dart_to_coamps' update the HDF5 files.
# 'trans_dart_to_coamps' uses COAMPS write routines that expect a DTG in
# the file name, the dart call provides the base filename of 'CoampsUpdate'
# Yes, I know the DTG is hardcoded here. My bad.

set instance = 1
while ( $instance <= ${ENSEMBLE_SIZE} )

   set  INFILE = `printf dart_%03d_output.nc         $instance`
   set OUTFILE = `printf coamps_%03d_2013011000.hdf5 $instance`

   cp ${TEMPLATE} ${OUTFILE}

   ln -sf ${INFILE}  dart_vector.nc
   ln -sf ${OUTFILE} CoampsUpdate_2013011000.hdf5

   ./trans_dart_to_coamps || exit 5

   @ instance ++

end

#===============================================================================

cat << README >! README.txt
The original file that was perturbed: ${ORIGIN}/${TEMPLATE}
README

#===============================================================================
echo
echo "${JOBNAME} ($JOBID) ended   at "`date`
echo
#===============================================================================

# optional cleanup block

# rm *.nc dart_log.* output_list*txt filter perfect_model_obs 
# rm trans_dart_to_coamps trans_coamps_to_dart dart.kstart
# rm coamps.hdf5 CoampsUpdate_2013011000.hdf5
# rm ${TEMPLATE}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
