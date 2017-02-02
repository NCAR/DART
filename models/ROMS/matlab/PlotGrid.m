function bob = PlotGrid(fname)
%% PlotGrid  plots the different ROMS grids and returns all the coordinate arrays.
%
% fname = '../data/wc13_grd.nc';
% romsmeta = PlotGrid(fname);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

ncid = netcdf.open(fname,'NC_NOWRITE');

% units are meters

x_rho = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'x_rho'));
x_u   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'x_u'));
x_v   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'x_v'));
x_psi = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'x_psi'));

y_rho = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'y_rho'));
y_u   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'y_u'));
y_v   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'y_v'));
y_psi = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'y_psi'));

% units are degrees

lon_rho = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lon_rho'));
lon_u   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lon_u'));
lon_v   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lon_v'));
lon_psi = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lon_psi'));

lat_rho = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lat_rho'));
lat_u   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lat_u'));
lat_v   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lat_v'));
lat_psi = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'lat_psi'));

% mask_psi:option_0 = "land" ;
% mask_psi:option_1 = "water" ;

mask_rho = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'mask_rho'));
mask_u   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'mask_u'));
mask_v   = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'mask_v'));
mask_psi = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'mask_psi'));

netcdf.close(ncid)

bob = struct( 'x_rho',    x_rho,    'x_u',    x_u,    'x_v',    x_v,    'x_psi', x_psi, ...
              'y_rho',    y_rho,    'y_u',    y_u,    'y_v',    y_v,    'y_psi', y_psi, ...
            'lon_rho',  lon_rho,  'lon_u',  lon_u,  'lon_v',  lon_v,  'lon_psi', lon_psi, ...
            'lat_rho',  lat_rho,  'lat_u',  lat_u,  'lat_v',  lat_v,  'lat_psi', lat_psi, ...
           'mask_rho', mask_rho, 'mask_u', mask_u, 'mask_v', mask_v, 'mask_psi', mask_psi);

figure(1); clf; orient tall

subplot(2,2,1)
plot(x_rho(:),y_rho(:),'x')
title('rho')

subplot(2,2,2)
plot(x_u(:),y_u(:),'x')
title('u')

subplot(2,2,3)
plot(x_v(:),y_v(:),'x')
title('v')

subplot(2,2,4)
plot(x_psi(:),y_psi(:),'x')
title('psi')

figure(2); clf; orient tall

subplot(2,2,1)
plot(lon_rho(:),lat_rho(:),'x')
title('rho')

subplot(2,2,2)
plot(lon_u(:),lat_u(:),'x')
title('u')

subplot(2,2,3)
plot(lon_v(:),lat_v(:),'x')
title('v')

subplot(2,2,4)
plot(lon_psi(:),lat_psi(:),'x')
title('psi')

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
