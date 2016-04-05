function h4 = plot_grid(fname)
%% plot_grid ... plots the ULAT,ULON and TLAT,TLON variables from a netcdf file.
% 
% fname = '../data/DATA.ionos2.nc';
% plot_grid(fname)
%

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

cg_lon     = nc_varget(fname,'cg_lon');
cg_lat     = nc_varget(fname,'cg_lat');
[cgx, cgy] = meshgrid(cg_lon, cg_lat);

ig_lon     = nc_varget(fname,'ig_lon');
ig_lat     = nc_varget(fname,'ig_lat');
[igx, igy] = meshgrid(ig_lon, ig_lat);

geo_lon = nc_varget(fname,'geo_lon');
geo_lat = nc_varget(fname,'geo_lat');

figure(1); clf; orient landscape

h1 = plot(cgx,cgy,'ko'); hold on;
h2 = plot(igx,igy,'rx');
h3 = plot(geo_lon,geo_lat,'gd');
title({'geographic is black o','magnetic is red x', 'mag2geo is green d'})
xlabel('longitude')
ylabel('latitude')
set(gca,'FontSize',20)

hold off;

%% Now plot just the first ROW
figure(2); clf; orient landscape

size(igx)
size(geo_lon)

h2 = plot(igx(:,1),igy(:,1),'rx'); hold on;
h3 = plot(geo_lon(:,1),geo_lat(:,1),'gd');
title({'ig_lon,ig_lat (magnetic) is red x', ...
       'geo_lon,geo_lat (mag2geo) is green d', ...
       'FIRST COLUMN ONLY','element (1,1) is black'},'Interpreter','none')
xlabel('longitude')
ylabel('latitude')
h4 = plot(igx(1,1),igy(1,1),'kx','LineWidth',2,'MarkerFaceColor','k');
h5 = plot(geo_lon(1,1),geo_lat(1,1),'kd','LineWidth',2,'MarkerFaceColor','k');
set(gca,'FontSize',20)

hold off;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

