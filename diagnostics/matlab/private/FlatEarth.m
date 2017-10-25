function FlatEarth(plotstruct)
%%
% plotstruct.axis    [xmin xmax ymin ymax zmin zmax] - all finite
% plotstruct.clim    [minimum data value, maximum data value]
% plotstruct.Ztype   a DART integer describing the vertical coordinate system

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------------
% GET THE ELEVATION DATA AND SET UP THE ASSOCIATED COORDINATE DATA
%---------------------------------------------------------------------------

topo = load('topo');     % GET Matlab-native [180x360] ELEVATION DATASET
lats = -89.5:89.5;       % CREATE LAT ARRAY FOR TOPO MATRIX
lons = 0.5:359.5;        % CREATE LON ARRAY FOR TOPO MATRIX
nlon = length(lons);
nlat = length(lats);

%%--------------------------------------------------------------------------
% IF WE NEED TO SWAP HEMISPHERES, DO SO NOW.
% If we didn't explicitly tell it, make a guess.
%---------------------------------------------------------------------------

axis(plotstruct.axis);
ax = plotstruct.axis;

if (ax(1) < 0)
   lons = lons - 180.0;
   topo.topo = [ topo.topo(:,nlon/2+1:nlon) topo.topo(:,1:nlon/2) ];
end

%%--------------------------------------------------------------------------
% We need to determine the geographic subset of the elevation matrix.
%---------------------------------------------------------------------------

lon_ind1 = find(ax(1) <= lons, 1);
lon_ind2 = find(ax(2) <= lons, 1);
lat_ind1 = find(ax(3) <= lats, 1);
lat_ind2 = find(ax(4) <= lats, 1);

if (isempty(lon_ind1)), lon_ind1 = 1;    end;
if (isempty(lon_ind2)), lon_ind2 = nlon; end;
if (isempty(lat_ind1)), lat_ind1 = 1;    end;
if (isempty(lat_ind2)), lat_ind2 = nlat; end;

elev = topo.topo(lat_ind1:lat_ind2,lon_ind1:lon_ind2);
x    = lons(lon_ind1:lon_ind2);
y    = lats(lat_ind1:lat_ind2);

%%--------------------------------------------------------------------------
% Augment the colormap and the CLim so that the lowest color index can be
% forced to a light gray without compromising the data range.

[clim, ~] = add_gray_color(plotstruct.clim);
cmin = double(clim(1)); % corresponds to the bogus data value

%%--------------------------------------------------------------------------
% Contour the "subset" - and give the whole thing an appropriate zlevel
% so the continents are either at the top of the plot (for depth), or
% the bottom (for just about everything else).
%@>todo FIXME ... observations AT the zlevel can be hidden under the landmasses
% when rotated to perfect overhead view.

orgholdstate = ishold; hold on;

set(gca,'Zdir','normal')
zlevel = ax(5);

if     (plotstruct.Ztype == -2)  %    zunits = 'undefined';
elseif (plotstruct.Ztype == -1)  %    zunits = 'surface';
elseif (plotstruct.Ztype ==  1)  %    zunits = 'level';
elseif (plotstruct.Ztype ==  2)  %    zunits = 'pressure';
    zlevel = ax(6);               % maximum pressure
    set(gca,'Zdir','reverse')
elseif (plotstruct.Ztype ==  3)  %    zunits = 'height';
end

% This translates to the right zlevel (no x,y translations).
% To set the color of the contours to the lowest color in the colormap
% add in the value that maps to the lowest color - the bogus data value.

if ( verLessThan('matlab','8.4') )
    [~,h] = contourf(x,y,elev+cmin,[cmin cmin],'k-');
    set(gca,'CLim',clim)
    t1    = hgtransform('Parent',gca);
    set(h,'Parent',t1);
    m     = makehgtform('translate',[0 0 zlevel]);
    set(t1,'Matrix',m)
else
    htrans = hgtransform;
    htrans.Matrix = makehgtform('translate',[0 0 zlevel]);
    [~,h] = contourf(x,y,elev+cmin,[cmin cmin],'k-','Parent',htrans);
    set(gca,'CLim',clim)
end
    
if (orgholdstate == 0), hold off; end;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
