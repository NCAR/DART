function plot_grid(fname)
%% plot_grid ... plots the ULAT,ULON and TLAT,TLON variables from a netcdf file.
% 
% fname = 'mpas_analysis.nc';
% plot_grid(fname)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

lat = nc_varget(fname,'latCell') * 180/pi;;
lon = nc_varget(fname,'lonCell') * 180/pi;;

inds = find(lon <= 0);
lon(inds) = lon(inds) + 360.0;

figure(1); clf; orient landscape
   ho = plot(lon,lat,'kx')
   axis([0 360 -90 90])
   title('MPAS cell centers')
   worldmap;
   axis image

figure(2); clf; orient landscape
   dt = DelaunayTri(lon,lat);
   voronoi(dt);
   axis([0 360 -90 90])
   title('MPAS cell centers and Voronoi diagram.')
   worldmap;
   axis image

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
