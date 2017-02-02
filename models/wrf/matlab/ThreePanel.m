%% ThreePanel - plots three graphics on one figure - uses Hui's data files.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clf; orient landscape; wysiwyg
clear all;

posmat = [0.250 0.5 0.500 0.5;
          0.075 0.1 0.400 0.4;
          0.525 0.1 0.400 0.4];
fontsize = 16;

subplot('position', posmat(1,:),'FontSize',fontsize)
fname = 'wrfinput_d01';

xlat      = nc_varget(fname,'XLAT');
xlon      = nc_varget(fname,'XLONG');
truelat1  = nc_attget(fname, nc_global, 'TRUELAT1');
truelat2  = nc_attget(fname, nc_global, 'TRUELAT2');
cen_lat   = nc_attget(fname, nc_global, 'CEN_LAT');
stand_lon = nc_attget(fname, nc_global, 'STAND_LON');

minlon = min(xlon(:)); maxlon = max(xlon(:));
minlat = min(xlat(:)); maxlat = max(xlat(:));

ax = axesm('lambert','Origin',[0 stand_lon 0],    ...
       'MapParallels',[truelat1 truelat2],        ...
       'MapLatLimit',[minlat maxlat], ...
       'MapLonLimit',[minlon maxlon] );

% this sets the limits for the entire projected domain
% [xlim ylim] = mfwdtran([xlat(1,1) xlat(end,end)],[xlon(1,1) xlon(end,end)]);
% set(gca,'xlim',xlim);
% set(gca,'ylim',ylim);

% This restricts the view to a subset of the domain
%set(gca,'xlim',[-0.1162 0.3084]);
%set(gca,'ylim',[ 0.1385 0.4089]);

% mlabel(16);
% plabel(-70);
% mlabel('on')

load ernesto2418.txt;
load ctl25.txt;
load gps25.txt;

lat_obs = ernesto2418(:, 3);
lon_obs = ernesto2418(:, 4);

lat_ctl = ctl25(:, 5);
lon_ctl = ctl25(:, 6);

lat_gps = gps25(:, 5);
lon_gps = gps25(:, 6);

plotm( lat_obs,lon_obs,'rx--', 'markersize', 8.0, 'LineWidth', 2.0);
plotm( lat_ctl,lon_ctl,'bx--', 'markersize', 8.0, 'LineWidth', 2.0);
plotm( lat_gps,lon_gps,'gv--', 'markersize', 8.0, 'LineWidth', 2.0);
legend('OBS', 'CTL', 'GPS'); legend boxoff

%coloring continent, slow
%geoshow(ax, xlat, xlon, 'DisplayType', 'texturemap', 'FaceColor', [1. 1. 1.])
bbox = zeros(2,2);
bbox(:) = [minlon minlat; maxlon maxlat];
landareas = shaperead('landareas', 'UseGeoCoords', true, ...
             'BoundingBox', bbox);
h = geoshow([landareas.Lat],[landareas.Lon],'Color','black');
% geoshow('landareas.shp', 'FaceColor', 'black');

%plotm( lat_non,lon_non,'g.--', 'markersize', 8.0);

%----------------------------------------------------------------------
% Track Error axis
%----------------------------------------------------------------------

subplot('position', posmat(2,:),'FontSize',fontsize)

lat_obs = ernesto2418(:, 3)*3.1415/180.0;
lon_obs = ernesto2418(:, 4)*3.1415/180.0;
int_obs = ernesto2418(:, 5);

lat_ctl = ctl25(:, 5)*3.1415/180.0;
lon_ctl = ctl25(:, 6)*3.1415/180.0;
int_ctl = ctl25(:, 7);

lat_gps = gps25(:, 5)*3.1415/180.0;
lon_gps = gps25(:, 6)*3.1415/180.0;
int_gps = gps25(:, 7);

%--------
lon_dist = lon_ctl - lon_obs;
for i= 1: size(lon_ctl);
ctl_dist(i) = acos(sin(lat_ctl(i)) * sin(lat_obs(i)) + cos(lat_ctl(i)) * cos(lat_obs(i)) * cos(lon_dist(i)))*6400.0;
end

lon_dist = lon_gps - lon_obs;
for i= 1: size(lon_ctl);
gps(i) = acos(sin(lat_gps(i)) * sin(lat_obs(i)) + cos(lat_gps(i)) * cos(lat_obs(i)) * cos(lon_dist(i)))*6400.0;
end

dd = ernesto2418(:,2)/2400.0  + ernesto2418(:,1);

plot( dd, ctl_dist,'bx--', dd, gps,'gv--',  'LineWidth', 2.0)
axis([25 28 0 150])
grid
xlabel('Day of August 2006)');
ylabel('Track error (km)');

legend('CTRL', 'GPS',  'Location', 'NorthWest'); legend boxoff



%----------------------------------------------------------------------
% Track Error axis
%----------------------------------------------------------------------

dd = ernesto2418(:,2)/2400.0  + ernesto2418(:,1);

subplot('position', posmat(3,:),'FontSize',fontsize)

plot( dd, int_obs,'rx-', dd, int_ctl,'bx--', dd, int_gps,'gv--', ...
      'LineWidth', 2.0);
%  change here
axis([25 28 990 1010])
set(gca,'YAxisLocation','right')
grid
xlabel('Day of August 2006');
ylabel('Surface Pressure Minimum (hPa)');
%legend('OBS', 'CTRL', 'GPS', 'Location', 'East')
legend('OBS', 'CTL', 'GPS','Location','SouthWest'); legend boxoff


print -dpng track25.png

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
