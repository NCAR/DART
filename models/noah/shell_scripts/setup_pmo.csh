#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This is an example script for how to stage the files in CENTRALDIR
# in preparation for a 'perfect model' or OSSE.
#
#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
#==============================================================================

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

#==============================================================================
# Stage all the required files in CENTRALDIR
#==============================================================================

set CENTRALDIR = `pwd`

set NOAHDIR = /Users/thoar/svn/DART/devel/models/noah/src/hrldas-v3.3
set DARTDIR = /Users/thoar/svn/DART/devel/models/noah

${COPY} ${NOAHDIR}/Run/Noah_hrldas_beta            .  || exit 1
${COPY} ${NOAHDIR}/Run/SOILPARM.TBL                .  || exit 1
${COPY} ${NOAHDIR}/Run/VEGPARM.TBL                 .  || exit 1
${COPY} ${NOAHDIR}/Run/GENPARM.TBL                 .  || exit 1
${COPY} ${NOAHDIR}/Run/URBPARM.TBL                 .  || exit 1

${COPY} ${DARTDIR}/templates/namelist.hrldas.template  namelist.hrldas || exit 2
${COPY} ${DARTDIR}/templates/wrfinput.template         wrfinput        || exit 2

${COPY} ${DARTDIR}/work/obs_seq.in                 .  || exit 3
${COPY} ${DARTDIR}/work/input.nml                  .  || exit 3
${COPY} ${DARTDIR}/work/perfect_model_obs          .  || exit 3
${COPY} ${DARTDIR}/work/dart_to_noah               .  || exit 3
${COPY} ${DARTDIR}/work/noah_to_dart               .  || exit 3
${COPY} ${DARTDIR}/shell_scripts/run_pmo.csh       .  || exit 3
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh .  || exit 3

# Need a single noah restart file to be used as THE TRUTH.
# This can come from anywhere. Link the file to that expected by DART:
# input.nml:model_nml noah_netcdf_filename = 'restart.nc'
# The assimilate.csh scripts wants an ensemble member node, for
# a perfect_model experiment - this must be 0001.

${COPY} ${DARTDIR}/ensemble_source/RESTART.2009010200_DOMAIN1.0001.nc  .  || exit 4
ln -sv RESTART.2009010200_DOMAIN1.0001.nc  restart.nc      || exit 5
ln -sv RESTART.2009010200_DOMAIN1.0001.nc  restart.0001.nc || exit 5


echo
echo "CENTRALDIR is ${CENTRALDIR}"
echo "Configure     ${CENTRALDIR}/input.nml"
echo "Configure     ${CENTRALDIR}/namelist.hrldas"
echo "Configure     ${CENTRALDIR}/wrfinput"
echo "execute       ${CENTRALDIR}/run_pmo.csh"
echo

exit 0


