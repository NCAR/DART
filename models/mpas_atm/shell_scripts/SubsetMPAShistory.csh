#!/bin/csh
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
# xtime                    coordinate variable
# latCell                  coordinate variable
# lonCell                  coordinate variable
# zgrid                    coordinate variable
# cellsOnVertex            coordinate variable
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

set MPASFILE = /gpfs/ptmp/syha/temp/mpas_init.nc
set MPASFILE = mpas_output.2010-10-23_03:00:00.nc

set COMPACT = mpas_analysis.nc

ncks -O -a -v xtime,latCell,lonCell,zgrid,cellsOnVertex,theta,uReconstructZonal,uReconstructMeridional,w,qv,qc,qr,qi,qs,qg,surface_pressure,edgeNormalVectors,nEdgesOnCell,edgesOnCell,u,dzs,xCell,yCell,zCell ${MPASFILE} ${COMPACT}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

