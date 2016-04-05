%% Visualize_Grid plots the grid, land mask, and tile fraction from a JULES static_init.nc

%% DART software - Copyright 2004 - 2015 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Visualize a 3x5 (nlon-X-nlat) grid file.

fname = '../work/static_input.nc';
ncid  = netcdf.open(fname,'NOWRITE');

lonid  = netcdf.inqVarID(ncid,'longitude');
lonmat = netcdf.getVar(ncid,lonid,'double');
latid  = netcdf.inqVarID(ncid,'latitude');
latmat = netcdf.getVar(ncid,latid,'double');

maskid   = netcdf.inqVarID(ncid,'land_mask');
landmask = netcdf.getVar(ncid,maskid,'double');
% landmask(1,5) = 0.5; for testing purposes.

tileid   = netcdf.inqVarID(ncid,'tile_frac');
tilefrac = netcdf.getVar(ncid,tileid,'double');

netcdf.close(ncid)

% This simple case can be represented by a single vector of latitudes and 
% another of longitudes.

lon = lonmat(:,1) ;
lat = latmat(1,:)';

lonleft  = lon(1) - (lon(2) - lon(1))/2;
lonright = lon(3) + (lon(3) - lon(2))/2;
lonedges = sort([lonleft mean(lon(1:2)) mean(lon(2:3)) lonright]);

lattop   = lat(1) + (lat(1) - lat(2))/2;
latbot   = lat(5) - (lat(4) - lat(5))/2;
latedges = sort([lattop mean(lat(1:2)) mean(lat(2:3)) mean(lat(3:4)) mean(lat(4:5)) latbot]);

[x,y] = meshgrid(lonedges,latedges);

clf; orient landscape
imagesc(lon,lat,landmask');
set(gca,'FontSize',14)
hold on;
plot(x,y,'x');
plot(lonmat,latmat,'o');

ax = axis;
ax = ax + [-0.001 0.001 -0.001 0.001];
axis(ax)

set(gca,'YTick',latedges)
set(gca,'XTick',lonedges)
grid on

[numx, numy] = size(lonmat);

for iy = 1:numy
for ix = 1:numx

   if (ix == 2 & iy == 3)
      mycolor = 'w';
   else
      mycolor = 'k';
   end

   str = sprintf('%d,%d',ix,iy);
   h = text(lonmat(ix,iy), latmat(ix,iy), str);
   set(h,'VerticalAlignment'  ,'bottom', 'Color', mycolor, ...
         'HorizontalAlignment','center','FontSize',18);

   str = sprintf('numtiles = %d',sum(tilefrac(ix,iy,:)> 0.0));
   h = text(lonmat(ix,iy), latmat(ix,iy), str);
   set(h,'VerticalAlignment'  ,'top', 'Color', mycolor, ...
         'HorizontalAlignment','center','FontSize',12);

end
end

% represent better precision on the X and Y axes

Xlabels = get(gca,'XTickLabels');
Ylabels = get(gca,'YTickLabels');

for ix = 1:numx+1
   Xlabels{ix} = sprintf('%f',lonedges(ix)+360.0);
%  Xlabels{ix} = sprintf('%f',lonedges(ix));
end

for iy = 1:numy+1
   Ylabels{iy} = sprintf('%f',latedges(iy));
end

set(gca, 'XTickLabels', Xlabels)
set(gca, 'YTickLabels', Ylabels)
set(gca, 'YDir', 'normal')
axis tight

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

