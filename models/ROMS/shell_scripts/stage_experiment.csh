#!/bin/tcsh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
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
# |-- input.nml
# |-- varinfo.dat
# |-- s4dvar.in
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
set EXPERIMENTDIR = /glade/scratch/${USER}/romstest2

if (-e ${EXPERIMENTDIR} ) then
   echo "ERROR: ${EXPERIMENTDIR} already exists."
   echo "Intentionally leaving it alone."
   echo "Must provide a new directory name for the experiement."
   exit 1
endif

#--------------------------------------------------------------------------
# stage everything in the experiment directory and then cat out some 
# instructions and put those same instructions in the experiment directory.
#--------------------------------------------------------------------------

mkdir -p ${EXPERIMENTDIR}
cd ${EXPERIMENTDIR}

rsync -Cavz ${ROMSDIR}/WC13/Data/      ${EXPERIMENTDIR}/Data/ || exit 1
rsync -Cavz ${ROMSDIR}/WC13/Ensemble/  ${EXPERIMENTDIR}/      || exit 1

\cp ${ROMSDIR}/External/varinfo.dat    ${EXPERIMENTDIR}/.     || exit 1
\cp ${ROMSDIR}/WC13/timtest2/oceanM    ${EXPERIMENTDIR}/.     || exit 1

\cp ${DARTDIR}/observations/ROMS/work/convert_roms_obs ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/work/input.nml              ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/work/filter                 ${EXPERIMENTDIR}/. || exit 2
\cp ${DARTDIR}/models/ROMS/shell_scripts/run_1day.csh  ${EXPERIMENTDIR}/. || exit 2

cat << ENDOFFILE >! README_assumptions.txt

#--------------------------------------------------------------------------
The EXPERIMENT directory is ${EXPERIMENTDIR}
The DART parts came from    ${DARTDIR}
The ROMS parts came from    ${ROMSDIR}
#--------------------------------------------------------------------------
The following changes must be made to ocean.in:

    NtileI*NtileJ must equal the task count given when ROMS was compiled

    NtileI == 4            1 yellowstone node has 16 procs, use them all
    NtileJ == 4            

    NTIMES == 48           \
        DT == 1800.0d0     ---  these two will advance ROMS for 1 day
      NRST == 48           write out a restart file after 1 day

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
and run the ${EXPERIMENTDIR}/run_1day.csh script.

ENDOFFILE

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size` || exit 3
set ensemble_size = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`         || exit 3

set member = 1
while ( ${member} <= ${ensemble_size} )

   set dirname = `printf instance_%04d $member`
   mkdir -p $dirname
   cd $dirname

   set UNIQUEFILENAME = `printf wc13_ini_%04d.nc $member`

   \cp ../ocean.in             .             || exit 4
   \mv ../${UNIQUEFILENAME}    .             || exit 4
   ln -sf ${UNIQUEFILENAME} roms_input.nc

   cd ..

   @ member++
end

#==========================================================================
# Then we run DART on the ensemble of new states
#==========================================================================

cat README_assumptions.txt


# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

