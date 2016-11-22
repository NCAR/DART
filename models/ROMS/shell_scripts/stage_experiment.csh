#!/bin/tcsh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
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
# |-- wc13_obs.nc
# |-- Data
# |   |-- adsen.cdl
# |   |-- <snip> whatever your configuration requires <snip>
# |   |-- wc13_grd.nc
# |   |-- wc13_ini.nc
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

set DARTDIR = /glade/p/work/${USER}/DART/rma_trunk
set ROMSDIR = /glade/p/work/${USER}/roms/test
set EXPERIMENTDIR = /glade/scratch/${USER}/roms_cycling_test
set SUBSTITUTE = /glade/p/work/${USER}/roms/trunk/ROMS/Bin/substitute

if (-e ${EXPERIMENTDIR} ) then
   echo "ERROR: ${EXPERIMENTDIR} already exists."
   echo "Intentionally leaving it alone."
   echo "Must provide a new directory name for the experiement."
   exit 1
endif

#--------------------------------------------------------------------------
# stage everything in the experiment directory
#--------------------------------------------------------------------------

mkdir -p ${EXPERIMENTDIR}
cd ${EXPERIMENTDIR}

rsync -Cavz ${ROMSDIR}/WC13/Data/      ${EXPERIMENTDIR}/Data/ || exit 1
rsync -Cavz ${ROMSDIR}/WC13/Ensemble/  ${EXPERIMENTDIR}/      || exit 1

\cp ${ROMSDIR}/External/varinfo.dat    ${EXPERIMENTDIR}/.     || exit 1
\cp ${ROMSDIR}/WC13/timtest2/oceanM    ${EXPERIMENTDIR}/.     || exit 1
\cp ${ROMSDIR}/WC13/timtest2/oceanG    ${EXPERIMENTDIR}/.

\cp ${DARTDIR}/models/ROMS/shell_scripts/run_multiple_jobs.csh  ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/shell_scripts/cycle.csh.template  ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/work/input.nml.template           ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/work/filter                       ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/observations/ROMS/work/convert_roms_obs       ${EXPERIMENTDIR}/. || exit 2

echo "no preexisting inflation files" >! ${EXPERIMENTDIR}/roms_inflation_cookie

# the observation file had some quirks

set ROMS_OBS = Data/obs_wc13_merged_2013+07d_pp1_depthinmeters_dt5400_physonly.nc

ncatted -a    units,survey_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
        -a calendar,survey_time,c,c,'gregorian' \
        -a    units,obs_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
        -a calendar,obs_time,c,c,'gregorian' \
        -a flag_meanings,obs_provenance,m,c,'gridded_AVISO_sea_level_anomaly_(zeta) gridded_Aquarius_SSS_(salinity) XBT_from_Met_Office_(temperature) CTD_from_Met_Office_(temperature) CTD_from_Met_Office_(salinity) ARGO_floats_(temperature) ARGO_floats_(salinity) glider_UCSD_(temperature) glider_UCSD_(salinity) blended_satellite_SST_(temperature)' \
        ${ROMS_OBS}

ncrename -a obs_type@long,long_name ${ROMS_OBS}

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

set ENSEMBLE_SIZE = 32
set ROMS_STDIN = ocean.in
set ROMS_DAPAR = s4dvar.in
set ROMS_DAI = roms_dai.nc
set ROMS_MOD = roms_mod_obs.nc
set ROMS_RST = roms_rst.nc
set ROMS_EXE = oceanM

# Set DART and ROMS (static) input values.
# There are some replacement strings left in the templates
# that will need to be replaced after each assimilation cycle.
# Things like DSTART and the ININAME ...

$SUBSTITUTE  ocean.in.template  MyVARNAME   ../varinfo.dat
$SUBSTITUTE  ocean.in.template  MyNtileI    4
$SUBSTITUTE  ocean.in.template  MyNtileJ    4
$SUBSTITUTE  ocean.in.template  MyNTIMES    48
$SUBSTITUTE  ocean.in.template  MyDT        3600.0d0
$SUBSTITUTE  ocean.in.template  MyNRST      24
$SUBSTITUTE  ocean.in.template  MyTIME_REF  19000101.0d0
$SUBSTITUTE  ocean.in.template  MyDAINAME   $ROMS_DAI
$SUBSTITUTE  ocean.in.template  MyRSTNAME   $ROMS_RST
$SUBSTITUTE  ocean.in.template  MyAPARNAM   $ROMS_DAPAR

$SUBSTITUTE  cycle.csh.template  MySUBSTITUTE          $SUBSTITUTE
$SUBSTITUTE  cycle.csh.template  EXPERIMENT_DIRECTORY  $EXPERIMENTDIR
$SUBSTITUTE  cycle.csh.template  MyROMS_EXE            $ROMS_EXE
$SUBSTITUTE  cycle.csh.template  MyROMS_STDIN          $ROMS_STDIN
$SUBSTITUTE  cycle.csh.template  MyOBSname             $ROMS_OBS
$SUBSTITUTE  cycle.csh.template  MyMODname             $ROMS_MOD
$SUBSTITUTE  cycle.csh.template  MyRSTNAME             $ROMS_RST
$SUBSTITUTE  cycle.csh.template  MyDAINAME             $ROMS_DAI

$SUBSTITUTE  run_multiple_jobs.csh  EXPERIMENT_DIRECTORY  $EXPERIMENTDIR

$SUBSTITUTE  s4dvar.in.template  MyMODname   $ROMS_MOD

$SUBSTITUTE  input.nml.template  Myens_size  $ENSEMBLE_SIZE
$SUBSTITUTE  input.nml.template  MyDAINAME   $ROMS_DAI

\cp ${DARTDIR}/system_simulation/final_full_precomputed_tables/final_full.$ENSEMBLE_SIZE .
\mv cycle.csh.template cycle.csh
chmod u+x cycle.csh
chmod u+x run_multiple_jobs.csh

set member = 1
while ( ${member} <= ${ENSEMBLE_SIZE} )

   set dirname = `printf instance_%04d $member`
   mkdir -p $dirname
   cd $dirname

   set ROMS_INI = `printf wc13_ini_%04d_2013_01_01.nc $member`

   \cp ../ocean.in.template    $ROMS_STDIN || exit 4
   \cp ../s4dvar.in.template   $ROMS_DAPAR || exit 4
   \mv ../$ROMS_INI            .           || exit 4

   # Set DSTART for the current ensemble.
   # NOTE ... bc can handle the 'long' integers that happen when the
   # reference time is 1900-01-01, the shell divide cannot.

   set OCEAN_TIME = `ncdump -v ocean_time ${ROMS_INI} | grep "ocean_time =" | tail -1`
   set TIME_SEC = `echo $OCEAN_TIME | grep -oE '[[:digit:]]+'`
   set DSTART = `echo "scale=6; $TIME_SEC / 86400.0" | bc `

   # Set ROMS standard input parameters needed in template scripts.

   $SUBSTITUTE $ROMS_STDIN MyDSTART   $DSTART
   $SUBSTITUTE $ROMS_STDIN MyININAME  $ROMS_INI

   $SUBSTITUTE $ROMS_DAPAR MyOBSname  ../$ROMS_OBS

   cd ..

   @ member++
end

# remove any leftover ensemble members
\rm wc13_ini_*.nc

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
   perturb_from_single_instance = .false.     (USUALLY!)
   use_restart_list             = .true.
   restart_list_file            = 'restart_files.txt'
   overwrite_state_input        = .true.
   obs_sequence_in_name         = "obs_seq.out"
   obs_sequence_out_name        = "obs_seq.final"
   ens_size                     = <YOUR ACTUAL ENSEMBLE SIZE>

&convert_roms_obs_nml
   roms_mod_obs_filelist        = 'precomputed_files.txt'
   dart_output_obs_file         = 'obs_seq.out'
   append_to_existing           = .false.
   use_precomputed_values       = .true.
   locations_in_IJK             = .false.

&model_nml
   output_state_vector          = .false.
   assimilation_period_days     = <something bigger than NTIMES*DT>
   assimilation_period_seconds  = 0
   vert_localization_coord      = 3

#--------------------------------------------------------------------------

After these files are confirmed, it should be possible to
run the ${EXPERIMENTDIR}/cycle.csh script.

ENDOFFILE

cat README_assumptions.txt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

