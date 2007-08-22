function pinfo = GetBgridInfo(pinfo_in,fname,routine);
% GetBgridInfo   prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetBgridInfo(pinfo_in,fname,routine);
%
% pinfo_in  Name of existing pinfo struct, e.g. output from CheckModelCompatibility
% fname     Name of the DART netcdf file
% routine   name of subsequent plot routine.

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

pinfo      = pinfo_in;
ft         = netcdf(fname);
model      = ft.model(:);
close(ft)

if strcmp(lower(model),'fms_bgrid') ~= 1
   error('Not so fast, this is not a bgrid model.')
end

copy   = getnc(fname,'copy');
times  = getnc(fname,'time');
levels = getnc(fname,'lev');
TmpI   = getnc(fname,'TmpI');    % temperature/pressure grid longitude
TmpJ   = getnc(fname,'TmpJ');    % temperature/pressure grid latitude
VelI   = getnc(fname,'VelI');    % velocity grid longitude
VelJ   = getnc(fname,'VelJ');    % velocity grid latitude

prognostic_vars  = {'ps','t','u','v'};

switch lower(deblank(routine))

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar           = GetVar(prognostic_vars);  % Determine prognostic variable
      [level, lvlind] = GetLevel(pgvar,levels);   % Determine level and index
      [lat  , latind] = GetLatitude( pgvar,TmpJ,VelJ);
      [lon  , lonind] = GetLongitude(pgvar,TmpI,VelI);

      pinfo = setfield(pinfo, 'model', model);
      pinfo = setfield(pinfo, 'fname', fname);
      pinfo = setfield(pinfo, 'var', pgvar);
      pinfo = setfield(pinfo, 'level', level);
      pinfo = setfield(pinfo, 'levelindex', lvlind);
      pinfo = setfield(pinfo, 'longitude', lon);
      pinfo = setfield(pinfo, 'lonindex', lonind);
      pinfo = setfield(pinfo, 'latitude', lat);
      pinfo = setfield(pinfo, 'latindex',latind);

   case 'plotcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(prognostic_vars);
      [base_time, base_tmeind] = GetTime(     base_var,times);
      [base_lvl,  base_lvlind] = GetLevel(    base_var,levels);
      [base_lat,  base_latind] = GetLatitude( base_var,TmpJ,VelJ);
      [base_lon,  base_lonind] = GetLongitude(base_var,TmpI,VelI);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvl);

      pinfo = setfield(pinfo, 'model', model);
      pinfo = setfield(pinfo, 'fname', fname);
      pinfo = setfield(pinfo, 'base_var', base_var);
      pinfo = setfield(pinfo, 'comp_var', comp_var);
      pinfo = setfield(pinfo, 'base_time', base_time);
      pinfo = setfield(pinfo, 'base_tmeind', base_tmeind);
      pinfo = setfield(pinfo, 'base_lvl', base_lvl);
      pinfo = setfield(pinfo, 'base_lvlind', base_lvlind);
      pinfo = setfield(pinfo, 'base_lat', base_lat);
      pinfo = setfield(pinfo, 'base_latind', base_latind);
      pinfo = setfield(pinfo, 'base_lon', base_lon);
      pinfo = setfield(pinfo, 'base_lonind', base_lonind);
      pinfo = setfield(pinfo, 'comp_lvl', comp_lvl);
      pinfo = setfield(pinfo, 'comp_lvlind', comp_lvlind);

   case 'plotvarvarcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(prognostic_vars);
      [base_time, base_tmeind] = GetTime(     base_var,times);
      [base_lvl , base_lvlind] = GetLevel(    base_var,levels);
      [base_lat , base_latind] = GetLatitude( base_var,TmpJ,VelJ);
      [base_lon , base_lonind] = GetLongitude(base_var,TmpI,VelI);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvl);
      [comp_lat, comp_latind] = GetLatitude( comp_var,TmpJ,VelJ, base_lat);
      [comp_lon, comp_lonind] = GetLongitude(comp_var,TmpI,VelI, base_lon);

      pinfo = setfield(pinfo, 'model', model);
      pinfo = setfield(pinfo, 'fname', fname);
      pinfo = setfield(pinfo, 'base_var', base_var);
      pinfo = setfield(pinfo, 'comp_var', comp_var);
      pinfo = setfield(pinfo, 'base_time', base_time);
      pinfo = setfield(pinfo, 'base_tmeind', base_tmeind);
      pinfo = setfield(pinfo, 'base_lvl', base_lvl);
      pinfo = setfield(pinfo, 'base_lvlind', base_lvlind);
      pinfo = setfield(pinfo, 'base_lat', base_lat);
      pinfo = setfield(pinfo, 'base_latind', base_latind);
      pinfo = setfield(pinfo, 'base_lon', base_lon);
      pinfo = setfield(pinfo, 'base_lonind', base_lonind);
      pinfo = setfield(pinfo, 'comp_lvl', comp_lvl);
      pinfo = setfield(pinfo, 'comp_lvlind', comp_lvlind);
      pinfo = setfield(pinfo, 'comp_lat', comp_lat);
      pinfo = setfield(pinfo, 'comp_latind', comp_latind);
      pinfo = setfield(pinfo, 'comp_lon', comp_lon);
      pinfo = setfield(pinfo, 'comp_lonind', comp_lonind);

   case 'plotsawtooth'

       pgvar          = GetVar(prognostic_vars);
      [level, lvlind] = GetLevel(pgvar,levels);   % Determine level and index
      [lat  , latind] = GetLatitude( pgvar,TmpJ,VelJ);
      [lon  , lonind] = GetLongitude(pgvar,TmpI,VelI);
      %[copy , lonind] = GetCopies(pgvar,copy);
      copyindices     = SetCopyID(fname);
      copy            = length(copyindices);

      pinfo = setfield(pinfo, 'model'         , model);
      pinfo = setfield(pinfo, 'var_names'     , pgvar);
      pinfo = setfield(pinfo, 'truth_file'    , []);
      %pinfo = setfield(pinfo, 'prior_file'    , pinfo.prior_file);
      %pinfo = setfield(pinfo, 'posterior_file', pinfo.posterior_file);
      pinfo = setfield(pinfo, 'level'         , level);
      pinfo = setfield(pinfo, 'levelindex'    , lvlind);
      pinfo = setfield(pinfo, 'latitude'      , lat);
      pinfo = setfield(pinfo, 'latindex'      , latind);
      pinfo = setfield(pinfo, 'longitude'     , lon);
      pinfo = setfield(pinfo, 'lonindex'      , lonind);
      pinfo = setfield(pinfo, 'copies'        , copy);
      pinfo = setfield(pinfo, 'copyindices'   , copyindices);

   case 'plotphasespace'

      disp('Getting information for the ''X'' variable.')
       var1                   = GetVar(prognostic_vars);
      [var1_lvl, var1_lvlind] = GetLevel(    var1, levels);
      [var1_lat, var1_latind] = GetLatitude( var1, TmpJ, VelJ);
      [var1_lon, var1_lonind] = GetLongitude(var1, TmpI, VelI);

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVar(prognostic_vars,        var1    );
      [var2_lvl, var2_lvlind] = GetLevel(    var2, levels,     var1_lvl);
      [var2_lat, var2_latind] = GetLatitude( var2, TmpJ, VelJ, var1_lat);
      [var2_lon, var2_lonind] = GetLongitude(var2, TmpI, VelI, var1_lon);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVar(prognostic_vars,        var1    );
      [var3_lvl, var3_lvlind] = GetLevel(    var3, levels,     var1_lvl);
      [var3_lat, var3_latind] = GetLatitude( var3, TmpJ, VelJ, var1_lat);
      [var3_lon, var3_lonind] = GetLongitude(var3, TmpI, VelI, var1_lon);

      % query for ensemble member
      s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
      if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''  ','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

      pinfo = setfield(pinfo, 'model',model);
      pinfo = setfield(pinfo, 'fname' ,fname);
      pinfo = setfield(pinfo, 'var1name' ,var1);
      pinfo = setfield(pinfo, 'var2name' ,var2);
      pinfo = setfield(pinfo, 'var3name' ,var3);
      pinfo = setfield(pinfo, 'var1_lvl' , var1_lvl);
      pinfo = setfield(pinfo, 'var1_lvlind', var1_lvlind);
      pinfo = setfield(pinfo, 'var1_lat' , var1_lat);
      pinfo = setfield(pinfo, 'var1_latind', var1_latind);
      pinfo = setfield(pinfo, 'var1_lon' , var1_lon);
      pinfo = setfield(pinfo, 'var1_lonind', var1_lonind);
      pinfo = setfield(pinfo, 'var2_lvl' , var2_lvl);
      pinfo = setfield(pinfo, 'var2_lvlind', var2_lvlind);
      pinfo = setfield(pinfo, 'var2_lat' , var2_lat);
      pinfo = setfield(pinfo, 'var2_latind', var2_latind);
      pinfo = setfield(pinfo, 'var2_lon' , var2_lon);
      pinfo = setfield(pinfo, 'var2_lonind', var2_lonind);
      pinfo = setfield(pinfo, 'var3_lvl' , var3_lvl);
      pinfo = setfield(pinfo, 'var3_lvlind', var3_lvlind);
      pinfo = setfield(pinfo, 'var3_lat' , var3_lat);
      pinfo = setfield(pinfo, 'var3_latind', var3_latind);
      pinfo = setfield(pinfo, 'var3_lon' , var3_lon);
      pinfo = setfield(pinfo, 'var3_lonind', var3_lonind);
      pinfo = setfield(pinfo, 'ens_mem'  , ens_mem);
      pinfo = setfield(pinfo, 'ltype',ltype);

   otherwise

end

function pgvar = GetVar(prognostic_vars, defvar)
%----------------------------------------------------------------------
if (nargin == 2), pgvar = defvar; else pgvar = 'ps'; end

str = sprintf(' %s ',prognostic_vars{1});
for i = 2:length(prognostic_vars),
   str = sprintf(' %s %s ',str,prognostic_vars{i});
end
disp(sprintf('Default variable is ''%s'', if this is OK, <cr>;',pgvar))  
disp(sprintf('If not, please enter one of: %s',str))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), pgvar = deblank(varstring); end 



function [time, timeind] = GetTime(pgvar, times, deftime)
%----------------------------------------------------------------------
if (nargin == 3), time = deftime; else time = mean(times); end

disp(sprintf('Default time is %f, if this is OK, <cr>;',time))
disp(sprintf('If not, enter a time between %.4f and %.4f, we use the closest.', ...
                         min(times),max(times)))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), time  = str2num(varstring); end 

d       = abs(time - times);    % crude distance
ind     = find(min(d) == d);  % multiple minima possible 
timeind = ind(1);             % use the first one
time    = times(timeind);


function [level, lvlind] = GetLevel(pgvar, levels, deflevel)
%----------------------------------------------------------------------
% level and lvlind will not be equal for all models, (and probably
% shouldn't for the bgrid ... but for future expansion ...
if (nargin == 3), level = deflevel; else level = 1; end

if strcmp(pgvar,'ps') ==1 
   disp('''ps'' only has one level, using it.')
   level  = 1;
   lvlind = 1;
else
   disp(sprintf('Default level is  %d, if this is OK, <cr>;',level))
   disp(sprintf('If not, enter a level between %d and %d, inclusive ...', ...
                         min(levels),max(levels)))
   varstring = input('we''ll use the closest (no syntax required)\n','s');

   if ~isempty(varstring), level = str2num(varstring); end 

   d      = abs(level - levels);  % crude distance
   ind    = find(min(d) == d);    % multiple minima possible 
   lvlind = ind(1);               % use the first one
   level  = levels(lvlind);
end



function [lon, lonind] = GetLongitude(pgvar, TmpI, VelI, deflon)
%----------------------------------------------------------------------
if (nargin == 4), lon = deflon; else lon = 255.0; end

switch lower(pgvar)
   case {'ps','t'}
     lons = TmpI;
   otherwise
     lons = VelI;
end

disp(sprintf('Default longitude is %f, if this is OK, <cr>;',lon))  
disp(sprintf('If not, enter a longitude between %.2f and %.2f, we use the closest.', ...
                         min(lons),max(lons)))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lon  = str2num(varstring); end 

d      = abs(lon - lons);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible 
lonind = ind(1);             % use the first one
lon    = lons(lonind);



function [lat, latind] = GetLatitude(pgvar, TmpJ, VelJ, deflat)
%----------------------------------------------------------------------
if (nargin == 4), lat = deflat; else lat = 40.0; end

switch lower(pgvar)
   case {'ps','t'}
     lats = TmpJ;
   otherwise
     lats = VelJ;
end

disp(sprintf('Default latitude is %f, if this is OK, <cr>;',lat))
disp(sprintf('If not, enter a latitude between %.2f and %.2f, we use the closest.', ...
                         min(lats),max(lats)))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lat = str2num(varstring); end 

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
