#!/bin/csh
#
# A set of examples on how to simply tweak existing netCDF files.
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

rm temp.nc lonlatlev.nc allgood.nc

set fname = /home/thoar/svn/DART/UCOAM/models/GCOM/data/InputGeo.nc

# remove the unwanted (state) variables from the geometry file

ncks -h -x -v u,v,w,time $fname temp.nc

# rename some variables and dimensions

ncrename -h -v depth,lev -v x,lon -v y,lat temp.nc
ncrename -h -d  lev,nzp1 -d  lat,nyp1 -d  lon,nxp1  \
            -d slev,nz   -d slat,ny   -d slon,nx  temp.nc

# manipulate the attributes

ncatted -h -a long_name,ulon,o,c,longitude temp.nc
ncatted -h -a long_name,vlon,o,c,longitude temp.nc
ncatted -h -a long_name,wlon,o,c,longitude temp.nc
ncatted -h -a long_name,lon,o,c,longitude temp.nc
ncatted -h -a axis,ulon,o,c,X temp.nc
ncatted -h -a axis,vlon,o,c,X temp.nc
ncatted -h -a axis,wlon,o,c,X temp.nc
ncatted -h -a axis,lon,o,c,X temp.nc
ncatted -h -a units,ulon,o,c,degrees_east temp.nc
ncatted -h -a units,vlon,o,c,degrees_east temp.nc
ncatted -h -a units,wlon,o,c,degrees_east temp.nc
ncatted -h -a units,lon,o,c,degrees_east temp.nc

ncatted -h -a long_name,ulat,o,c,latitude temp.nc
ncatted -h -a long_name,vlat,o,c,latitude temp.nc
ncatted -h -a long_name,wlat,o,c,latitude temp.nc
ncatted -h -a long_name,lat,o,c,latitude temp.nc
ncatted -h -a axis,ulat,o,c,Y temp.nc
ncatted -h -a axis,vlat,o,c,Y temp.nc
ncatted -h -a axis,wlat,o,c,Y temp.nc
ncatted -h -a axis,lat,o,c,Y temp.nc
ncatted -h -a units,ulat,o,c,degrees_north temp.nc
ncatted -h -a units,vlat,o,c,degrees_north temp.nc
ncatted -h -a units,wlat,o,c,degrees_north temp.nc
ncatted -h -a units,lat,o,c,degrees_north temp.nc

ncatted -h -a long_name,ulev,o,c,depth temp.nc
ncatted -h -a long_name,vlev,o,c,depth temp.nc
ncatted -h -a long_name,wlev,o,c,depth temp.nc
ncatted -h -a long_name,lev,o,c,depth temp.nc
ncatted -h -a axis,ulev,o,c,Z temp.nc
ncatted -h -a axis,vlev,o,c,Z temp.nc
ncatted -h -a axis,wlev,o,c,Z temp.nc
ncatted -h -a axis,lev,o,c,Z temp.nc
ncatted -h -a units,ulev,o,c,meters temp.nc
ncatted -h -a units,vlev,o,c,meters temp.nc
ncatted -h -a units,wlev,o,c,meters temp.nc
ncatted -h -a units,lev,o,c,meters temp.nc
ncatted -h -a positive,ulev,o,c,up temp.nc
ncatted -h -a positive,vlev,o,c,up temp.nc
ncatted -h -a positive,wlev,o,c,up temp.nc
ncatted -h -a positive,lev,o,c,up temp.nc

# subset correct variables to a temp file

ncks -h -x -v lon,lev,lat temp.nc allgood.nc

# subset lon(nzp1,nyp1,nxp1) to proper shape
# i.e.   lon(nz,ny,nx) to proper shape
# by skipping dimension 0

ncks -a -h -d nzp1,1, -dnyp1,1, -dnxp1,1, -v lon,lev,lat temp.nc lonlatlev.nc

# rename the dimensions to nz,ny,nx

ncrename -h -d nzp1,nz -d  nyp1,ny -d nxp1,nx lonlatlev.nc

# paste the two files into one. Finis.

ncks -a -h -A lonlatlev.nc allgood.nc

#=======================================================================

# This part removes the geometry variables from the restart file

set fname = gcom_restart_original.nc

rm temp.nc

ncks -h -v time,u,v,w,p,T,S $fname temp.nc
ncatted -h -a coordinates,u,o,c,'ulon ulat ulev time' temp.nc
ncatted -h -a coordinates,v,o,c,'vlon vlat vlev time' temp.nc
ncatted -h -a coordinates,w,o,c,'wlon wlat wlev time' temp.nc
ncatted -h -a coordinates,p,o,c,'lon lat lev time' temp.nc
ncatted -h -a coordinates,S,o,c,'lon lat lev time' temp.nc
ncatted -h -a coordinates,T,o,c,'lon lat lev time' temp.nc

ncatted -h -a long_name,time,o,c,time temp.nc
ncatted -h -a calendar,time,o,c,gregorian temp.nc
ncatted -h -a units,time,o,c,'days since 1986-01-01 00:00:00' temp.nc

ncatted -h -a units,p,o,c,bar temp.nc
ncatted -h -a units,S,o,c,psu temp.nc
ncatted -h -a units,T,o,c,'degrees C' temp.nc
ncatted -h -a units,u,o,c,'m/s' temp.nc
ncatted -h -a units,v,o,c,'m/s' temp.nc
ncatted -h -a units,w,o,c,'m/s' temp.nc

ncrename -h -d nx,nxp1 -d ny,nyp1 -d nz,nzp1 temp.nc
ncrename -h -d sx,nx   -d sy,ny   -d sz,nz   temp.nc

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

