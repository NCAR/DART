%% gloc.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

nc=netcdf('wrfinput_d01');

xlat=nc{'XLAT'}(1,:,:);
xlon=nc{'XLONG'}(1,:,:);
truelat1 = nc.TRUELAT1(1);
truelat2 = nc.TRUELAT2(1);
%cen_lat = nc.CEN_LAT(1);
%cen_lon = nc.CEN_LON(1);
stand_lon = nc.STAND_LON;

load coast;

ax=axesm('lambert','origin', [0 stand_lon 0], 'mapparallels', [truelat1 truelat2], ...
      'MapLatLimit',[min(min(xlat(:))) max(max(xlat(:)))], ...
      'MapLonLimit', [min(min(xlon(:))) max(max(xlon(:)))],...
       'fontcolor', [.1 .1 .1]);
%      'fontcolor', [.1 .1 .1], ...
%      'mlinelocation',15,'plinelocation',15,...
%      'meridianlabel','on','mlabellocation',15,...
%      'parallellabel','on','plabellocation',15,...
%      'fontcolor', [.75 .75 .75]);

%gridm; plotm(lat,long,'color',[.75 .75 .75]); tightmap
gridm; plotm(lat,long,'color',[1 1 1]); tightmap
%mlabel; plabel;
mlabel(20); plabel(-130);
%plabel('fontcolor',[0.75 0.75 0.75]);

[xlim ylim]=mfwdtran([xlat(1,1) xlat(end,end)],[xlon(1,1) xlon(end,end)]);
set(gca,'xlim',[min(xlim(:)) max(xlim(:))]);
set(gca,'ylim',[min(ylim(:)) max(ylim(:))]);

%coloring continent, slow
geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [1. 1. 1.])
%geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [.96 .91 .76])
%-------------------------------------------------------------


load loc_2007jun08.dat;
lat_pts = loc_2007jun08(:, 1);
lon_pts = loc_2007jun08(:, 2);

% subplot('position', [0.1,0.2,0.7,0.7]);
 plotm(lat_pts,lon_pts,'k.',  'markersize', 2.0);

%plotm(lat_pts,lon_pts,'b.', 'markersize', 8.8);

print -dps gloc_08jun2007.ps
print -dpng gloc_08jun2007.png

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
