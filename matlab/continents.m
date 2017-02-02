function h = continents(ifill,lonorg,level)
%% CONTINENTS  overlays the continents on an existing figure as a F(lon,lat).
%           The continents can be filled or "hollow".
%
%           This is not appropriate for use with a real MAP (i.e. anything
%           other than a simple lon vs. lat plot)
%
% USAGE: continents([fill],[lonorigin],[level])
%
% fill == 'hollow' 		=> outline of landmasses (DEFAULT) 
% fill == 'solid' 		=> solid black landmasses
% fill == 'light' 		=> solid grey landmasses
% lonorigin == 'greenwich'	=> Greenwich Meridian  [   0,360]
% lonorigin == 'dateline'	=> Greenwich Meridian  [-180,180]
% level == 123                  => Z location to plot the 2d sheet
%
% The order of the arguments is immaterial. They are all optional.
%
% This basically gets 1degree elevation data and then contours it
% around zero ... It tries to be intelligent about whether you are 
% plotting [-180,180] or [0,360], but you can override that with "lonorigin".
%
% EXAMPLE:
%
% load topo;				% GET [180x360] ELEVATION DATASET
% lats = [-89.5:89.5];			% CREATE LAT ARRAY FOR TOPO MATRIX
% lons = [0.5:359.5];			% CREATE LON ARRAY FOR TOPO MATRIX
% imagesc(lons,lats,topo);		% CREATE SOME PLOT w/ true x,y limits
% set(gca,'YDir','normal');		% CORRECT ORIENTATION
%
% continents;				% COASTLINES ONLY -- OR --
%
% continents('solid');			% SOLID LANDMASSES

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

filltype  = 'hollow';
lonorigin = 'greenwich';
zlevel = 0.0;

if (nargin >= 1)

   if (isa(ifill,'char'))
      switch lower(ifill)
         case 'light', 		filltype  = 'light';
         case 'hollow', 	filltype  = 'hollow';
         case 'solid',  	filltype  = 'solid';
         case 'greenwich',	lonorigin = 'greenwich';
         case 'dateline',  	lonorigin = 'dateline';
         otherwise, error('argument unknown')
      end
   else
      zlevel = ifill;
   end

end

if (nargin >= 2)

   if (isa(lonorg,'char'))
      switch lower(lonorg)
         case 'light', 		filltype  = 'light';
         case 'hollow', 	filltype  = 'hollow';
         case 'solid',  	filltype  = 'solid';
         case 'greenwich',	lonorigin = 'greenwich';
         case 'dateline',  	lonorigin = 'dateline';
         otherwise, error('second argument unknown')
      end
   else
      zlevel = lonorg;
   end

end
   
if (nargin >= 3)

   if (isa(level,'char'))
      switch lower(level)
         case 'light', 		filltype  = 'light';
         case 'hollow', 	filltype  = 'hollow';
         case 'solid',  	filltype  = 'solid';
         case 'greenwich',	lonorigin = 'greenwich';
         case 'dateline',  	lonorigin = 'dateline';
         otherwise, error('third argument unknown')
      end
   else
      zlevel = level;
   end

end
   
if (nargin >= 4)

   disp(' ')
   disp('AAAAAAAAHHHHH! Thump.')
   disp(' ')
   error('Wrong number of arguments.')

end

%---------------------------------------------------------------------------
% GET THE ELEVATION DATA AND SET UP THE ASSOCIATED COORDINATE DATA
%---------------------------------------------------------------------------

load topo;                 % GET Matlab-native [180x360] ELEVATION DATASET
lats = [-89.5:89.5];                    % CREATE LAT ARRAY FOR TOPO MATRIX
lons = [0.5:359.5];                     % CREATE LON ARRAY FOR TOPO MATRIX
nlon = length(lons);
nlat = length(lats);

topo = topo + zlevel;       % optional shift of vertical level 

%---------------------------------------------------------------------------
% IF WE NEED TO SWAP HEMISPHERES, DO SO NOW.
% If we didn't explicitly tell it, make an educated guess.
%---------------------------------------------------------------------------

ax   = axis;

if strcmp(lonorigin,'dateline')
   lons = lons - 180.0;
   topo = [ topo(:,nlon/2+1:nlon) topo(:,1:nlon/2) ];
elseif strcmp(lonorigin,'greenwich') 
%   disp('Greenwich on left')
%   do nothing.
elseif (ax(1) < -2) 
   lons = lons - 180.0;
   topo = [ topo(:,nlon/2+1:nlon) topo(:,1:nlon/2) ];
end

%---------------------------------------------------------------------------
% We need to determine the geographic subset of the elevation matrix.
%---------------------------------------------------------------------------

lon_ind1 = min(find(ax(1) <= lons));
lon_ind2 = min(find(ax(2) <= lons));
lat_ind1 = min(find(ax(3) <= lats));
lat_ind2 = min(find(ax(4) <= lats));
 
if (isempty(lon_ind1)) lon_ind1 = 1;    end;
if (isempty(lon_ind2)) lon_ind2 = nlon; end;
if (isempty(lat_ind1)) lat_ind1 = 1;    end;
if (isempty(lat_ind2)) lat_ind2 = nlat; end;
 
elev = topo(lat_ind1:lat_ind2,lon_ind1:lon_ind2);
x    = lons(lon_ind1:lon_ind2);
y    = lats(lat_ind1:lat_ind2);

%---------------------------------------------------------------------------
% Contour the "subset"
% There are differences between 6.5 and 7.0 that make changing the colors
% of the filled contours a real pain. Providing both solutions.
%---------------------------------------------------------------------------

orgholdstate = ishold;
hold on;

% outlines only.  contour3 works in both 2D and 3D plots.
if strcmp(filltype,'hollow')

   [c,h] = contour3(x,y,elev,[zlevel zlevel],'k-');

else

   % must use contourf first, then alter color and possibly level.
   if strcmp(filltype,'light')
      fcolor = [0.7 0.7 0.7];    % light grey
   else
      fcolor = [0.0 0.0 0.0];    % black
   end

   % compute the filled contour lines.
   [c,h] = contourf(x,y,elev,[zlevel zlevel],'k-');

   % give them the requested color
   if (length(h) == 1)    % matlab v 7.0 and above
      kids = get(h,'Children'); 
      set(kids,'FaceColor',fcolor);
   else                   % matlab 6.5 and below
      set(h   ,'FaceColor',fcolor);
   end

   % if not in the 2d plane, add z coordinates at the right level
   % (contourf is 2D only; no filled contour3() option exists.)
   if (zlevel ~= 0.0)  
      if (length(h) == 1)    % matlab v 7.0 and above
         kids = get(h,'Children'); 
         for i=1:length(kids)
           mz = get(kids(i), 'XData');
           z = zeros(length(mz), 1) + zlevel;
           set(kids(i), 'ZData', z);
         end
      else                   % matlab 6.5 and below
         for i=1:length(h)
            mz = get(h(i), 'XData');
            z = zeros(length(mz), 1) + zlevel;
            set(h(i), 'ZData', z);
         end
      end
   end

end

if (orgholdstate == 0) hold off; end;

axis(ax);				% MAINTAIN ORIGINAL LIMITS


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
