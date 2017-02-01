#!/bin/sh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script simply pares down an MPAS history file to a more 'tidy' size,
# at the expense of losing information we may need but currently do not 
# envision needing.
#
# DART needs the following coordinate variables:
# xtime                    time variable
# latCell                  cell center coordinate variable
# lonCell                  cell center coordinate variable
# zgrid                    height coordinate variable
# cellsOnVertex            define triangular mesh coordinate variable
# pressure_base            vertical coordinate variable (base)
# pressure_p               vertical coordinate variable (perturbation)
#
# theta                    potential state vector component
# uReconstructZonal        potential state vector component
# uReconstructMeridional   potential state vector component
# w                        potential state vector component
# qv                       potential state vector component
# qc                       potential state vector component
# qr                       potential state vector component
# qi                       potential state vector component
# qs                       potential state vector component
# qg                       potential state vector component
#
# surface_pressure         chosen because it is only (Time,nCells)
# edgeNormalVectors        needed by statevector_to_analysis_file()
# nEdgesOnCell             needed by statevector_to_analysis_file()
# u                        needed by statevector_to_analysis_file()
# dzs                      just to get the 'nSoilLevels' dimension

MPASFILE=/gpfs/ptmp/syha/temp/mpas_init.nc
#MPASFILE=mpas_output.2010-10-23_03:00:00.nc

COMPACT=mpas_analysis.nc

# i will happily put this back to csh if you can tell
# me how to specify a single token that's multiline.
# csh kept putting spaces in this variable between lines.
# i want a single 'word' here.
FLIST=xtime,latCell,lonCell,zgrid,cellsOnVertex,\
theta,pressure_base,pressure_p,exner,\
uReconstructZonal,uReconstructMeridional,\
rho,w,qv,qc,qr,qi,qs,qg,surface_pressure,edgeNormalVectors,\
nEdgesOnCell,edgesOnCell,u,dzs,xCell,yCell,zCell

echo copying over the following fields: $FLIST

ncks -O -a -v ${FLIST} ${MPASFILE} ${COMPACT}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

