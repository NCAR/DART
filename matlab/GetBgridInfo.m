function pinfo = GetBgridInfo(fname,routine);
% GetBgridInfo   prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetBgridInfo(fname,routine);
%
% fname     Name of the DART netcdf file
% routine   name of subsequent plot routine.

ft         = netcdf(fname);
model      = ft.model(:);
close(ft)

if strcmp(lower(model),'fms_bgrid') ~= 1
   error('Not so fast, this is not a bgrid model.')
end

copy   = getnc(fname,'copy');
times  = getnc(fname,'time');
levels = getnc(fname,'level');
TmpI   = getnc(fname,'TmpI');    % temperature/pressure grid longitude
TmpJ   = getnc(fname,'TmpJ');    % temperature/pressure grid latitude
VelI   = getnc(fname,'VelI');    % velocity grid longitude
VelJ   = getnc(fname,'VelJ');    % velocity grid latitude

prognostic_vars  = {'ps','t','u','v'};

switch lower(deblank(routine))

   case 'plotbins'

      pgvar           = GetVar(prognostic_vars);  % Determine prognostic variable
      [level, lvlind] = GetLevel(pgvar,levels);   % Determine level and index
      [lon  , lonind] = GetLongitude(pgvar,TmpI,VelI);
      [lat  , latind] = GetLatitude( pgvar,TmpJ,VelJ);

      pinfo = struct('model',model, ...
              'fname',fname, ...
              'var',pgvar, ...
              'level'    ,level, 'levelindex',lvlind, ...
              'longitude',lon  ,   'lonindex',lonind, ...
              'latitude' ,lat  ,   'latindex',latind);

   case 'plotcorrel'

      pgvar           = GetVar(prognostic_vars);  % Determine prognostic variable
      [level, lvlind] = GetLevel(pgvar,levels);   % Determine level and index
      [lon  , lonind] = GetLongitude(pgvar,TmpI,VelI);
      [lat  , latind] = GetLatitude( pgvar,TmpJ,VelJ);

      pinfo = struct('model',model, ...
              'fname',fname, ...
              'var',pgvar, ...
              'level'    ,level, 'levelindex',lvlind, ...
              'longitude',lon  ,   'lonindex',lonind, ...
              'latitude' ,lat  ,   'latindex',latind);

   otherwise

end

function pgvar = GetVar(prognostic_vars)
%----------------------------------------------------------------------
disp('Default variable is ''ps'', if this is OK, <cr>;')  
disp('If not, please enter one of following possible ')
prognostic_vars
varstring = input('(no syntax required)\n','s');
if isempty(varstring)
   pgvar = 'ps';
else
   pgvar = deblank(varstring);
end 


function [level, lvlind] = GetLevel(pgvar,levels)
%----------------------------------------------------------------------
% level and lvlind will not be equal for all models, (and probably
% shouldn't for the bgrid ... but for future expansion ...
level = NaN;

if strcmp(pgvar,'ps') ==1 
   disp('''ps'' only has one level, using it.')
   level  = 1;
   lvlind = 1;
else
   while all(level ~= levels)   % if you don't match any of the levels, try again
      disp('Default level is  1  , if this is OK, <cr>;')  
      disp(sprintf('If not, enter a level between %d and %d, inclusive.', ...
                         min(levels),max(levels)))
      varstring = input('(no syntax required)\n','s');
      if isempty(varstring)
         level = 1.0;
         lvlind = level;
      else
         level = str2num(varstring);
         lvlind = level;
      end 
   end % of the while loop.
end



function [lon, lonind] = GetLongitude(pgvar,TmpI,VelI)
%----------------------------------------------------------------------
lon = NaN;

switch lower(pgvar)
   case {'ps','t'}
     lons = TmpI;
   otherwise
     lons = VelI;
end

disp('Default longitude is 255 E (Boulder,CO) if this is OK, <cr>;')  
disp(sprintf('If not, enter a longitude between %.2f and %.2f, we use the closest.', ...
                         min(lons),max(lons)))
varstring = input('(no syntax required)\n','s');
if isempty(varstring)
   lon = 255.0;
else
   lon  = str2num(varstring);
end 
d      = abs(lon - lons);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible 
lonind = ind(1);             % use the first one
lon    = lons(lonind);



function [lat, latind] = GetLatitude(pgvar,TmpJ,VelJ)
%----------------------------------------------------------------------
lat = NaN;

switch lower(pgvar)
   case {'ps','t'}
     lats = TmpJ;
   otherwise
     lats = VelJ;
end

disp('Default latitude is 40 N (Boulder,CO) if this is OK, <cr>;')  
disp(sprintf('If not, enter a latitude between %.2f and %.2f, we use the closest.', ...
                         min(lats),max(lats)))
varstring = input('(no syntax required)\n','s');
if isempty(varstring)
   lat = 40.0;
else
   lat = str2num(varstring);
end 
d      = abs(lat - lats);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible 
latind = ind(1);             % use the first one
lat    = lats(latind);



function dist = arcdist(lat,lon,lat2,lon2)

% arcdist  great arc distance (km) between points on an earth-like sphere.
%
% function dist = arcdist(lat,lon,latvec,lonvec)
%
% lat,lon    MUST be scalars 
% lat1,lon1  lat1,lon1 can be any shape (but must match each other)
%
% enter: 
%       lat1,lon1 = lat, lon (in degrees)
%       lat2,lon2 = lat, lon (for loc 2)
% returns:
%       dist = [vector of] great arc distance between points
%
% Assumes a spherical earth and fractional lat/lon coords
%
% Example (1 degree at the equator):  
% dist = arcdist(0.0,0.0,0.0,1.0)	% 1 degree longitude
% dist = arcdist(0.0,0.0,1.0,0.0)       % 1 degree latitude
% dist = arcdist(60.0,0.0,60.0,1.0)     % 1 degree longitude at a high latitude
% dist = arcdist(40.0,-105.251,42.996,-78.84663)  % over the river and ...

r      = 6378.136;		% equatorial radius of earth in km

Dlat2 = size(lat2);
Dlon2 = size(lon2);

Lat2   = lat2(:);
Lon2   = lon2(:);

alpha  = abs(Lon2 - lon)*pi/180;
colat1 =    (90   - Lat2)*pi/180;
colat2 =    (90   - lat)*pi/180;

ang1   = cos(colat2).*cos(colat1);
ang2   = sin(colat2).*sin(colat1).*cos(alpha);

if ( prod(Dlat1) == 1 ) 
   dist   = acos(ang1 + ang2)*r;
else
   dist   = reshape(acos(ang1 + ang2)*r,Dlat2);
end
