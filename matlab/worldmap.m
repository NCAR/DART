function h = worldmap(ifill,lonorg,level)
% WORLDMAP  overlays the continents on an existing figure as a F(lon,lat).
%           The continents can be filled or "hollow".
%
%           This is not appropriate for use with a real MAP (i.e. anything
%           other than a simple lon vs. lat plot)
%
% USAGE: worldmap([fill],[lonorigin],[level])
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
% worldmap;				% COASTLINES ONLY -- OR --
%
% worldmap('solid');			% SOLID LANDMASSES
%

% tjh   1998.mar.26

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2006, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

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

if strcmp(filltype,'hollow')

   [c,h] = contour(x,y,elev,[zlevel zlevel],'k-');

elseif strcmp(filltype,'light')

   [c,h] = contourf(x,y,elev,[zlevel zlevel],'k-');

   if (length(h) == 1)    % do it the 'new' way
      kids = get(h,'Children'); 
      set(kids,'FaceColor',[0.7 0.7 0.7]); % SET COLOR TO lt gray
   else
      set(h   ,'FaceColor',[0.7 0.7 0.7]); % SET COLOR TO lt gray
   end

else

   [c,h] = contourf(x,y,elev,[zlevel zlevel],'k-');

   if (length(h) == 1)    % do it the 'new' way
      kids = get(h,'Children'); 
      set(kids,'FaceColor',[0.0 0.0 0.0]); % SET COLOR TO black
   else
      set(h   ,'FaceColor',[0.0 0.0 0.0]); % SET COLOR TO black
   end

end

if (orgholdstate == 0) hold off; end;

axis(ax);				% MAINTAIN ORIGINAL LIMITS

