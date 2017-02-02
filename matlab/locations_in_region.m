function inds = locations_in_region(locations, region)
%% locations_in_region returns the indices of the locations in the specified region.
%
% USAGE:
% inds = locations_in_region(locations, region);
%
% locations is N-by-3 array of locations - [lon; lat; lvl]
% region is an N-by-1 array of region boundaries:
%        [leftlon rightlon minlat maxlat [minz maxz]]
%        If leftlon is numerically larger than rightlon,
%        the region will span the Prime Meridian.
%
% EXAMPLE: simple case
% inds = locations_in_region(locations,[180 360 -90 90]); % Western Hemisphere
%
% EXAMPLE: longitude region spanning the Prime Meridian
% inds = locations_in_region(locations,[330 180 -90 90]); % Everything BUT the Americas
%
% EXAMPLE: vertical subsetting
% inds = locations_in_region(locations,[0 360 -90 90 0 1000]); % Everything  between 0,1000

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (length(region) == 6)
   zmin = min(region(5:6));
   zmax = max(region(5:6));
elseif (length(region) < 4)
   error('region must be an array of length 4 or 6')
end

ymin = min(region(3:4));
ymax = max(region(3:4));
xmin = mod(region(1),360.0);
xmax = mod(region(2),360.0);

lons = mod(locations(:,1),360.0);

% find the ones in the right latitudes & levels

latlogical = ((locations(:,2) >= ymin) & (locations(:,2) <= ymax));
lvllogical = ((locations(:,3) >= zmin) & (locations(:,3) <= zmax));

% find the ones in the right longitudes ... including wrapping

if ( xmin == xmax )
   lonlogical = ones(size(latlogical));   % everybody plays
else

   if ( xmin > xmax )
      xmax = xmax + 360.0;
      inds = find(lons < xmin);
      lons(inds) = lons(inds) + 360.0;
   end
   lonlogical = ((lons >= xmin) & (lons <= xmax));
end

% put it all together ...

inds = find( (lonlogical .* latlogical .* lvllogical) > 0 );

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
