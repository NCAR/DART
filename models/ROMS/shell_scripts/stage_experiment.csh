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

set DARTDIR = /glade/p/work/${USER}/DART/rma_roms
set ROMSDIR = /glade/p/work/${USER}/roms/test
set EXPERIMENTDIR = /glade/scratch/${USER}/romstest3
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

\cp ${DARTDIR}/models/ROMS/shell_scripts/cycle.csh.template  ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/work/input.nml.template     ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/work/filter                 ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/observations/ROMS/work/convert_roms_obs ${EXPERIMENTDIR}/. || exit 2

#--------------------------------------------------------------------------
# customize the user input 
# ROMS_STDIN_TMP and ROMS_DAPAR_TMP templates come from the Ensemble directory.
# The ocean.in.template and s4dvar.in.template files have hardcoded strings
# to be replaced using the SUBSTITUTE command. These templates are copied to
# each ROMS instance directory and modified with the values declared here.
#--------------------------------------------------------------------------

set ROMS_STDIN_TMP = ocean.in.template
set ROMS_DAPAR_TMP = s4dvar.in.template
set ROMS_STDIN = ocean.in
set ROMS_DAPAR = s4dvar.in
set ROMS_EXE = oceanM

set VARNAME = varinfo.dat
set NtileI = 4
set NtileJ = 4
set NCPUS = 16
set NTIMES = 48
set ROMS_DT = 3600.0d0
set NRST = 24
set TIME_REF = 19000101.0d0

set ROMS_DAI = roms_dai.nc
set ROMS_INI = roms_ini.nc
set ROMS_RST = roms_rst.nc
set ROMS_MOD = roms_mod.nc

# There are some values that DART needs so specify

set ENSEMBLE_SIZE = 3

\cp input.nml.template input.nml
\cp cycle.csh.template cycle.csh

#--------------------------------------------------------------------------
# put some instructions in the experiment directory and echo to screen
#--------------------------------------------------------------------------

cat << ENDOFFILE >! README_assumptions.txt

#--------------------------------------------------------------------------
The EXPERIMENT directory is ${EXPERIMENTDIR}
The DART parts came from    ${DARTDIR}
The ROMS parts came from    ${ROMSDIR}
#--------------------------------------------------------------------------
The following changes must be made to ocean.in:

    NtileI*NtileJ must equal the task count

    NtileI == 4            1 yellowstone node has 16 procs, use them all
    NtileJ == 4

    NTIMES == 48           \
        DT == 1800.0d0     ---  these two will advance ROMS for 1 day
      NRST == 48           write out a restart file after 1 day

    NTIMES == 24           \
        DT == 3600.0d0     ---  these two will advance ROMS for 1 day
      NRST == 24           write out a restart file after 1 day

   ININAME == roms_input.nc
   VARNAME == ../varinfo.dat
   APARNAM == ../s4dvar.in      required for writing VERIFICATION obs
   DAINAME == roms_dai.nc

#--------------------------------------------------------------------------
The following changes must be made to s4dvar.in

    OBSname == ../Data/roms_obs_in.nc     input observations
    MODname == roms_obs_mod.nc            ROMS-computed observations

#--------------------------------------------------------------------------
The following settings MUST be made to input.nml, there are many
more that you _may_ want to change or set.

&filter_nml:
   perturb_from_single_instance = .false.     (USUALLY!)
   use_restart_list             = .true.
   restart_list_file            = 'restart_files.txt'
   overwrite_state_input        = .true.
   obs_sequence_in_name         = "obs_seq.out"
   obs_sequence_out_name        = "obs_seq.final"
   ens_size                     = <YOUR ACTUAL ENSEMBLE SIZE>

&convert_roms_obs_nml
   roms_input_obs_file          = 'Data/wc13_obs.nc'
   roms_mod_obs_files           = ''
   roms_mod_obs_filelist        = 'precomputed_files.txt'
   dart_output_obs_file         = 'obs_seq.out'
   locations_in_IJK             = .false.
   add_random_noise             = .true.
   pert_amplitude               = 0.01
   append_to_existing           = .false.

&model_nml
   output_state_vector          = .false.
   assimilation_period_days     = <something bigger than NTIMES*DT>
   assimilation_period_seconds  = 0
   roms_filename                = 'roms_input.nc'
   vert_localization_coord      = 3

#--------------------------------------------------------------------------

After these changes are made, it should be possible to configure
and run the ${EXPERIMENTDIR}/cycle.csh script.

ENDOFFILE

set ROMS_OBS = Data/obs_wc13_merged_2013+07d_pp1_depthinmeters_dt5400_physonly.nc

ncatted -a    units,survey_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
        -a calendar,survey_time,c,c,'gregorian' \
        -a    units,obs_time,m,c,'days since 1900-01-01 00:00:00 GMT' \
        -a calendar,obs_time,c,c,'gregorian' \
        -a flag_meanings,obs_provenance,m,c,'gridded_AVISO_sea_level_anomaly_(zeta) gridded_Aquarius_SSS_(salinity) XBT_from_Met_Office_(temperature) CTD_from_Met_Office_(temperature) CTD_from_Met_Office_(salinity) ARGO_floats_(temperature) ARGO_floats_(salinity) glider_UCSD_(temperature) glider_UCSD_(salinity) blended_satellite_SST_(temperature)' \
        ${ROMS_OBS}

ncrename -a obs_type@long,long_name ${ROMS_OBS}

$SUBSTITUTE input.nml Myens_size            ${ENSEMBLE_SIZE}
$SUBSTITUTE input.nml MyMODname             ${ROMS_MOD}
$SUBSTITUTE input.nml MyDAINAME             ${ROMS_DAI}

$SUBSTITUTE cycle.csh EXPERIMENT_DIRECTORY  ${EXPERIMENTDIR}
$SUBSTITUTE cycle.csh MyROMS_EXE            ${ROMS_EXE}
$SUBSTITUTE cycle.csh MyROMS_STDIN          ${ROMS_STDIN}
$SUBSTITUTE cycle.csh MyDAINAME             ${ROMS_DAI}
$SUBSTITUTE cycle.csh MyININAME             ${ROMS_INI}
$SUBSTITUTE cycle.csh MyOBSname             ${ROMS_OBS}
$SUBSTITUTE cycle.csh MyMODname             ${ROMS_MOD}

set member = 1
while ( ${member} <= ${ENSEMBLE_SIZE} )

   set dirname = `printf instance_%04d $member`
   mkdir -p $dirname
   cd $dirname

   set ROMS_INI = `printf wc13_ini_%04d_2013_01_01.nc $member`

   \cp ../${ROMS_STDIN_TMP}    ${ROMS_STDIN} || exit 4
   \cp ../${ROMS_DAPAR_TMP}    ${ROMS_DAPAR} || exit 4
   \mv ../${ROMS_INI}          .             || exit 4

   # Set DSTART for the current ensemble.
   # NOTE ... bc can handle the 'long' integers that happen when the
   # reference time is 1900-01-01, the shell divide cannot.

   set OCEAN_TIME = `ncdump -v ocean_time ${ROMS_INI} | grep "ocean_time =" | tail -1`
   set TIME_SEC = `echo $OCEAN_TIME | grep -oE '[[:digit:]]+'`
   set DSTART = `echo "scale=6; $TIME_SEC / 86400.0" | bc `

   # Set ROMS standard input parameters needed in template scripts.

   $SUBSTITUTE $ROMS_STDIN MyVARNAME  ../$VARNAME
   $SUBSTITUTE $ROMS_STDIN MyNtileI   $NtileI
   $SUBSTITUTE $ROMS_STDIN MyNtileJ   $NtileJ
   $SUBSTITUTE $ROMS_STDIN MyNTIMES   $NTIMES
   $SUBSTITUTE $ROMS_STDIN MyDT       $ROMS_DT
   $SUBSTITUTE $ROMS_STDIN MyNRST     $NRST
   $SUBSTITUTE $ROMS_STDIN MyDSTART   $DSTART
   $SUBSTITUTE $ROMS_STDIN MyTIME_REF $TIME_REF
   $SUBSTITUTE $ROMS_STDIN MyININAME  $ROMS_INI
   $SUBSTITUTE $ROMS_STDIN MyDAINAME  $ROMS_DAI
   $SUBSTITUTE $ROMS_STDIN MyRSTNAME  $ROMS_RST
   $SUBSTITUTE $ROMS_STDIN MyAPARNAM  $ROMS_DAPAR

   $SUBSTITUTE $ROMS_DAPAR MyOBSname  $ROMS_OBS
   $SUBSTITUTE $ROMS_DAPAR MyMODname  $ROMS_MOD

   cd ..

   @ member++
end

# remove any leftover ensemble members
\rm wc13_ini_*.nc

#==========================================================================
# Then we run DART on the ensemble of new states
#==========================================================================

cat README_assumptions.txt

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

