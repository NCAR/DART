#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#==========================================================================
#
# This script prepares a directory for an assimilation experiment.
# The idea is basically that everything you need should be assembled
# by this script and that this should only be run ONCE per experiment.
# After everything is staged in the experiment directory, another script
# can be run to advance the model and perform the assimilation.
#
# The EXPERIMENTDIR should be structured something like:
#
# |-- oceanM
# |-- convert_roms_obs
# |-- filter
# |-- input.nml.template
# |-- varinfo.dat
# |-- s4dvar.in.template
# |-- ocean.in.template
# |-- wc12_obs.nc
# |-- Data
# |   |-- adsen.cdl
# |   |-- <snip> whatever your configuration requires <snip>
# |   |-- wc12_grd.nc
# |   |-- wc12_ini.nc
# |-- instance_0001
# |   |-- ocean.in
# |   `-- roms_input.nc
# |-- instance_...
# |   |-- ocean.in
# |   `-- roms_input.nc
# `-- instance_NNNN
#     |-- ocean.in
#     `-- roms_input.nc
#
# This script can be run interactively.
#==========================================================================

# DARTDIR declares the root location of the DART code tree.
# ROMSDIR declares the location of the ROMS 'test' code tree
#         because that's where I built my ROMS and I'm using the
#         default forcing/data files.

switch ("`hostname`")
   case eddy
      set DARTDIR = /home/${USER}/DART/rma_trunk
      set ROMSDIR = /home/${USER}/WC12_DART
      set  SRCDIR = /home/amm/ROMS/TRUNK_JAN17/ROMS
      set SUBSTITUTE = ${SRCDIR}/Bin/substitute
      set EXPERIMENTDIR = /home/${USER}/thoar_eddy1/roms_cycling_test
      breaksw
   case ch*
      set DARTDIR = /glade/p/work/${USER}/DART/rma_trunk
      set ROMSDIR = /glade/p/work/${USER}/roms/WC12
      set  SRCDIR = /glade/p/work/${USER}/roms/trunk/ROMS
      set SUBSTITUTE = ${SRCDIR}/Bin/substitute
      set EXPERIMENTDIR = /glade/scratch/${USER}/roms_cheyenne_test
      breaksw
   case ys*
      set DARTDIR = /glade/p/work/${USER}/DART/rma_trunk
      set ROMSDIR = /glade/p/work/${USER}/roms/WC12
      set  SRCDIR = /glade/p/work/${USER}/roms/trunk/ROMS
      set SUBSTITUTE = ${SRCDIR}/Bin/substitute
      set EXPERIMENTDIR = /glade/scratch/${USER}/roms_cycling_test
      breaksw
   default
      breaksw
endsw

if (-e ${EXPERIMENTDIR} ) then
   echo "ERROR: ${EXPERIMENTDIR} already exists."
   echo "Intentionally leaving it alone."
   echo "Must provide a new directory name for the experiement."
   exit 1
endif

#--------------------------------------------------------------------------
# stage everything in the experiment directory
#--------------------------------------------------------------------------

set ENSEMBLE_SIZE = 5
set ROMS_STDIN = ocean.in
set ROMS_DAPAR = s4dvar.in
set ROMS_DAI = roms_dai.nc
set ROMS_MOD = roms_mod_obs.nc
set ROMS_RST = roms_rst.nc
set ROMS_EXE = oceanM
set ROMS_OBS = Obs/obs_37623.nc
set ROMS_INIBASE = wc12_ini

mkdir -p ${EXPERIMENTDIR}
cd ${EXPERIMENTDIR}

rsync -Cavz ${ROMSDIR}/Obs/                 ${EXPERIMENTDIR}/Obs/  || exit 1

foreach FILE ( \
    ${ROMSDIR}/oceanM                                                     \
    ${ROMSDIR}/Ensemble/ocean.in.template                                 \
    ${ROMSDIR}/Ensemble/s4dvar.in.template                                \
    ${SRCDIR}/External/varinfo.dat                                        \
    ${DARTDIR}/observations/obs_converters/ROMS/work/convert_roms_obs     \
    ${DARTDIR}/models/ROMS/work/input.nml.template                        \
    ${DARTDIR}/models/ROMS/work/filter                                    \
    ${DARTDIR}/models/ROMS/shell_scripts/get_ocean_time.csh               \
    ${DARTDIR}/models/ROMS/shell_scripts/cycle.csh.template               \
    ${DARTDIR}/models/ROMS/shell_scripts/submit_multiple_cycles_lsf.csh   \
    ${DARTDIR}/models/ROMS/shell_scripts/advance_ensemble.csh.template    \
    ${DARTDIR}/models/ROMS/shell_scripts/run_filter.csh.template          \
    ${DARTDIR}/models/ROMS/shell_scripts/submit_multiple_jobs_slurm.csh   )
    \cp ${FILE} ${EXPERIMENTDIR}/. || exit 2
end

set SAMPDIR = assimilation_code/programs/system_simulation/work
set SAMPFILE = sampling_error_correction_table.Lanai.nc
\cp ${DARTDIR}/${SAMPDIR}/${SAMPFILE} sampling_error_correction_table.nc

echo "no preexisting inflation files" >! ${EXPERIMENTDIR}/roms_inflation_cookie

# THIS PARTICULAR observation file had some quirks.
# These changes should not be needed in general.
#
# set ROMS_OBS = Data/obs_wc13_merged_2013+07d_pp1_depthinmeters_dt5400_physonly.nc
#
# ncatted -a    units,survey_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
#         -a calendar,survey_time,c,c,'gregorian' \
#         -a    units,obs_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
#         -a calendar,obs_time,c,c,'gregorian' \
#         ${ROMS_OBS}

#--------------------------------------------------------------------------
# customize the user input templates with things that will remain constant
# througout the assimilation. We want to leave the *.template files with
# items that will change with each assimilation cycle.
#
# ocean.in.template and s4dvar.in.template come from the Ensemble directory.
# The ocean.in.template and s4dvar.in.template files have hardcoded strings
# to be replaced using the SUBSTITUTE command. These templates are copied to
# each ROMS instance directory and modified with the values declared here.
#--------------------------------------------------------------------------

foreach FILE ( ocean.in.template s4dvar.in.template )
   if ( ! -e $FILE ) then
      echo "ERROR: $FILE must exist."
      echo "ERROR: $FILE expected to come from the Ensemble directory."
      exit 1
   endif
end


# Set DART and ROMS (static) input values.
# There are some replacement strings left in the templates
# that will need to be replaced after each assimilation cycle.
# Things like DSTART and the ININAME ...

$SUBSTITUTE  ocean.in.template  MyVARNAME    ../varinfo.dat
$SUBSTITUTE  ocean.in.template  MyNtileI     4
$SUBSTITUTE  ocean.in.template  MyNtileJ     4
$SUBSTITUTE  ocean.in.template  MyNTIMES     96
$SUBSTITUTE  ocean.in.template  MyDT         900.0d0
$SUBSTITUTE  ocean.in.template  MyNRST       96
$SUBSTITUTE  ocean.in.template  MyTIME_REF   19000101.0d0
$SUBSTITUTE  ocean.in.template  MyDAINAME    $ROMS_DAI
$SUBSTITUTE  ocean.in.template  MyRSTNAME    $ROMS_RST
$SUBSTITUTE  ocean.in.template  MyAPARNAM    $ROMS_DAPAR
$SUBSTITUTE  ocean.in.template  MyROMS_STDIN $ROMS_STDIN

$SUBSTITUTE  s4dvar.in.template  MyMODname   $ROMS_MOD

$SUBSTITUTE  input.nml.template  Myens_size  $ENSEMBLE_SIZE
$SUBSTITUTE  input.nml.template  MyDAINAME   $ROMS_DAI

# submit_multiple_cycles_lsf.csh and cycle.csh go together.

$SUBSTITUTE  submit_multiple_cycles_lsf.csh  EXPERIMENT_DIRECTORY  $EXPERIMENTDIR
chmod u+x submit_multiple_cycles_lsf.csh

$SUBSTITUTE  cycle.csh.template  MySUBSTITUTE          $SUBSTITUTE
$SUBSTITUTE  cycle.csh.template  EXPERIMENT_DIRECTORY  $EXPERIMENTDIR
$SUBSTITUTE  cycle.csh.template  MyROMS_EXE            $ROMS_EXE
$SUBSTITUTE  cycle.csh.template  MyROMS_STDIN          $ROMS_STDIN
$SUBSTITUTE  cycle.csh.template  MyMODname             $ROMS_MOD
$SUBSTITUTE  cycle.csh.template  MyRSTNAME             $ROMS_RST
$SUBSTITUTE  cycle.csh.template  MyDAINAME             $ROMS_DAI
\mv cycle.csh.template cycle.csh
chmod u+x cycle.csh

# submit_multiple_jobs_slurm.csh, advance_ensemble.csh and run_filter.csh go together.

$SUBSTITUTE  submit_multiple_jobs_slurm.csh  EXPERIMENT_DIRECTORY  $EXPERIMENTDIR
chmod u+x submit_multiple_jobs_slurm.csh

$SUBSTITUTE  advance_ensemble.csh.template  Myens_size            $ENSEMBLE_SIZE
$SUBSTITUTE  advance_ensemble.csh.template  MySUBSTITUTE          $SUBSTITUTE
$SUBSTITUTE  advance_ensemble.csh.template  EXPERIMENT_DIRECTORY  $EXPERIMENTDIR
$SUBSTITUTE  advance_ensemble.csh.template  MyROMS_EXE            $ROMS_EXE
$SUBSTITUTE  advance_ensemble.csh.template  MyROMS_STDIN          $ROMS_STDIN
$SUBSTITUTE  advance_ensemble.csh.template  MyMODname             $ROMS_MOD
$SUBSTITUTE  advance_ensemble.csh.template  MyRSTNAME             $ROMS_RST
$SUBSTITUTE  advance_ensemble.csh.template  MyDAINAME             $ROMS_DAI
$SUBSTITUTE  advance_ensemble.csh.template  MyAPARNAM             $ROMS_DAPAR
$SUBSTITUTE  advance_ensemble.csh.template  MyINIBASE             $ROMS_INIBASE
\mv advance_ensemble.csh.template advance_ensemble.csh
chmod u+x advance_ensemble.csh

$SUBSTITUTE  run_filter.csh.template  MySUBSTITUTE          $SUBSTITUTE
$SUBSTITUTE  run_filter.csh.template  EXPERIMENT_DIRECTORY  $EXPERIMENTDIR
$SUBSTITUTE  run_filter.csh.template  MyDAINAME             $ROMS_DAI
$SUBSTITUTE  run_filter.csh.template  MyMODname             $ROMS_MOD
$SUBSTITUTE  run_filter.csh.template  MyROMS_STDIN          $ROMS_STDIN
\mv run_filter.csh.template run_filter.csh
chmod u+x run_filter.csh

set member = 1
while ( ${member} <= ${ENSEMBLE_SIZE} )

   set dirname = `printf instance_%04d $member`
   mkdir -p $dirname
   cd $dirname

   set ROMS_INI = `printf %s_%04d.nc $ROMS_INIBASE $member`

   \cp ../ocean.in.template           $ROMS_STDIN || exit 4
   \cp ../s4dvar.in.template          $ROMS_DAPAR || exit 4
   \cp ${ROMSDIR}/Ensemble/$ROMS_INI       .      || exit 4

   $SUBSTITUTE $ROMS_STDIN MyININAME  $ROMS_INI

   cd ..

   @ member++
end

#--------------------------------------------------------------------------
# put some instructions in the experiment directory and echo to screen
#--------------------------------------------------------------------------

cat << ENDOFFILE >! README_assumptions.txt

#--------------------------------------------------------------------------
The EXPERIMENT directory is ${EXPERIMENTDIR}
The DART parts came from    ${DARTDIR}
The ROMS parts came from    ${ROMSDIR}
#--------------------------------------------------------------------------

    NtileI*NtileJ must equal the task count

    NtileI == 4            1 yellowstone node has 16 procs, use them all
    NtileJ == 4

    NTIMES == 48           \
        DT == 1800.0d0     ---  these two will advance ROMS for 1 day
      NRST == 48           write out a restart file after 1 day

    NTIMES == 24           \
        DT == 3600.0d0     ---  these two will advance ROMS for 1 day
      NRST == 24           write out a restart file after 1 day

#--------------------------------------------------------------------------
The following settings MUST be made to input.nml.template for the logic
in cycle.csh to work correctly.  There are many more that you _may_ want
to change or set to impact the performance of the assimilation.

&filter_nml:
   input_state_file_list        = 'restart_files.txt'
   output_state_file_list       = 'restart_files.txt'
   obs_sequence_in_name         = 'obs_seq.out'
   obs_sequence_out_name        = 'obs_seq.final'
   ens_size                     = <YOUR ACTUAL ENSEMBLE SIZE>

&convert_roms_obs_nml
   roms_mod_obs_filelist        = 'precomputed_files.txt'
   dart_output_obs_file         = 'obs_seq.out'
   append_to_existing           = .false.
   use_precomputed_values       = .true.

&model_nml
   assimilation_period_days     = <something bigger than NTIMES*DT>
   assimilation_period_seconds  = 0
   vert_localization_coord      = 3

#--------------------------------------------------------------------------

After these files are confirmed, it should be possible to
run any of these scripts in ${EXPERIMENTDIR} :
   cycle.csh
   run_filter.csh
   submit_multiple_jobs_slurm.csh
.

ENDOFFILE

cat README_assumptions.txt

exit 0


