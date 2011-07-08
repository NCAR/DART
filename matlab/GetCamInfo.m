function pinfo = GetCamInfo(pstruct,routine);
%% GetCamInfo   prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetCamInfo(pstruct,routine);
%
% pstruct   structure containing the names of the truth_file and the diagn_file of the DART netcdf file
% routine   name of subsequent plot routine.

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (isfield(pstruct,'truth_file') && exist(pstruct.truth_file,'file') )
       fname = pstruct.truth_file;
elseif (isfield(pstruct,'diagn_file') && exist(pstruct.diagn_file,'file') )
       fname = pstruct.diagn_file;
elseif (isfield(pstruct,'prior_file') && exist(pstruct.prior_file,'file') )
       fname = pstruct.prior_file;
elseif (isfield(pstruct,'posterior_file') && exist(pstruct.posterior_file,'file') )
       fname = pstruct.posterior_file;
elseif (isfield(pstruct,'fname') && exist(pstruct.fname,'file') )
       fname = pstruct.fname;
end

model  = nc_attget(fname,nc_global,'model');

if strcmpi(model,'cam') ~= 1
   error('Not so fast, this is not a cam model.')
end

pinfo  = pstruct;

copy   = nc_varget(fname,'copy');
times  = nc_varget(fname,'time');
ilevel = nc_varget(fname,'ilev');    % interfaces
levels = nc_varget(fname, 'lev');    % midpoints
lon    = nc_varget(fname, 'lon');
lat    = nc_varget(fname, 'lat');

% A more robust way would be to use the netcdf low-level ops:
% bob = var(f);     % bob is a cell array of ncvars
% name(bob{1})       % is the variable name string
% bob{1}(:)          % is the value of the netcdf variable  (no
% offset/scale)
% have not yet figured out a way to only use non-coordinate variables.

prognostic_vars = {'PS','T','U','V','Q','CLDLIQ','CLDICE'};
num_vars = length(prognostic_vars);

switch lower(deblank(routine))

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar           = GetVar(prognostic_vars);  % Determine prognostic variable
      [level, lvlind] = GetLevel(pgvar,levels);   % Determine level and index
      [lat  , latind] = GetLatitude( pgvar,lat);
      [lon  , lonind] = GetLongitude(pgvar,lon);

      pinfo.model      = model;
      pinfo.var        = pgvar;
      pinfo.level      = level;
      pinfo.levelindex = lvlind;
      pinfo.longitude  = lon;
      pinfo.lonindex   = lonind;
      pinfo.latitude   = lat;
      pinfo.latindex   = latind;


   case 'plotcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(prognostic_vars);
      [base_time, base_tmeind] = GetTime(     base_var,times);
      [base_lvl,  base_lvlind] = GetLevel(    base_var,levels);
      [base_lat,  base_latind] = GetLatitude( base_var,lat);
      [base_lon,  base_lonind] = GetLongitude(base_var,lon);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvlind);

      pinfo.model       = model; 
      pinfo.base_var    = base_var;
      pinfo.comp_var    = comp_var;
      pinfo.base_time   = base_time;
      pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl    = base_lvl;
      pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat    = base_lat;
      pinfo.base_latind = base_latind;
      pinfo.base_lon    = base_lon;
      pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl    = comp_lvl;
      pinfo.comp_lvlind = comp_lvlind;


   case 'plotvarvarcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(prognostic_vars);
      [base_time, base_tmeind] = GetTime(     base_var,times);
      [base_lvl , base_lvlind] = GetLevel(    base_var,levels);
      [base_lat , base_latind] = GetLatitude( base_var,lat);
      [base_lon , base_lonind] = GetLongitude(base_var,lon);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvlind);
      [comp_lat, comp_latind] = GetLatitude( comp_var,lat,base_lat);
      [comp_lon, comp_lonind] = GetLongitude(comp_var,lon,base_lon);

      pinfo.model       = model;
      pinfo.base_var    = base_var;
      pinfo.comp_var    = comp_var;
      pinfo.base_time   = base_time;
      pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl    = base_lvl;
      pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat    = base_lat;
      pinfo.base_latind = base_latind;
      pinfo.base_lon    = base_lon;
      pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl    = comp_lvl;
      pinfo.comp_lvlind = comp_lvlind;
      pinfo.comp_lat    = comp_lat;
      pinfo.comp_latind = comp_latind;
      pinfo.comp_lon    = comp_lon;
      pinfo.comp_lonind = comp_lonind;


   case 'plotsawtooth'

       pgvar          = GetVar(prognostic_vars);
      [level, lvlind] = GetLevel(    pgvar,levels);
      [  lat, latind] = GetLatitude( pgvar,lat);
      [  lon, lonind] = GetLongitude(pgvar,lon);
 %    [  lon, lonind] = GetCopies(pgvar,xxx);

      pinfo.model          = model;
      pinfo.var_names      = pgvar;
      pinfo.truth_file     = [];
      pinfo.prior_file     = pstruct.prior_file;
      pinfo.posterior_file = pstruct.posterior_file;
      pinfo.level          = level;
      pinfo.levelindex     = lvlind;
      pinfo.latitude       = lat;
      pinfo.latindex       = latind;
      pinfo.longitude      = lon;
      pinfo.lonindex       = lonind;
      pinfo.copies         = 0;
      pinfo.copyindices    = [];

      if ( exist(pstruct.truth_file,'file') )
         pinfo.truth_file = pstruct.truth_file;
      end
       

   case 'plotphasespace'

      disp('Getting information for the ''X'' variable.')
       var1                   = GetVar(prognostic_vars);
      [var1_lvl, var1_lvlind] = GetLevel(    var1, levels);
      [var1_lat, var1_latind] = GetLatitude( var1, lat );
      [var1_lon, var1_lonind] = GetLongitude(var1, lon );

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVar(prognostic_vars,        var1    );
      [var2_lvl, var2_lvlind] = GetLevel(    var2, levels,     var1_lvlind);
      [var2_lat, var2_latind] = GetLatitude( var2, lat, var1_lat);
      [var2_lon, var2_lonind] = GetLongitude(var2, lon, var1_lon);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVar(prognostic_vars,        var1    );
      [var3_lvl, var3_lvlind] = GetLevel(    var3, levels,     var1_lvlind);
      [var3_lat, var3_latind] = GetLatitude( var3, lat, var1_lat);
      [var3_lon, var3_lonind] = GetLongitude(var3, lon, var1_lon);

      % query for ensemble member
      s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
      if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''  ','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

      pinfo.model       = model;
      pinfo.var1name    = var1;
      pinfo.var2name    = var2;
      pinfo.var3name    = var3;
      pinfo.var1_lvl    = var1_lvl;
      pinfo.var1_lvlind = var1_lvlind;
      pinfo.var1_lat    = var1_lat;
      pinfo.var1_latind = var1_latind;
      pinfo.var1_lon    = var1_lon;
      pinfo.var1_lonind = var1_lonind;
      pinfo.var2_lvl    = var2_lvl;
      pinfo.var2_lvlind = var2_lvlind;
      pinfo.var2_lat    = var2_lat;
      pinfo.var2_latind = var2_latind;
      pinfo.var2_lon    = var2_lon;
      pinfo.var2_lonind = var2_lonind;
      pinfo.var3_lvl    = var3_lvl;
      pinfo.var3_lvlind = var3_lvlind;
      pinfo.var3_lat    = var3_lat;
      pinfo.var3_latind = var3_latind;
      pinfo.var3_lon    = var3_lon;
      pinfo.var3_lonind = var3_lonind;
      pinfo.ens_mem     = ens_mem;
      pinfo.ltype       = ltype;

   otherwise

end


function pgvar = GetVar(prognostic_vars, defvar)
%----------------------------------------------------------------------
if (nargin == 2), pgvar = defvar; else pgvar = 'PS'; end

str = sprintf(' %s ',prognostic_vars{1});
for i = 2:length(prognostic_vars),
   str = sprintf(' %s %s ',str,prognostic_vars{i});
end
fprintf('Default variable is ''%s'', if this is OK, <cr>;\n',pgvar)  
fprintf('If not, please enter one of: %s\n',str)
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), pgvar = strtrim(varstring); end 
inds        = strfind(pgvar,',');
pgvar(inds) = '';


function [time, timeind] = GetTime(pgvar, times, deftime)
%----------------------------------------------------------------------
if (nargin == 3), time = deftime; else time = mean(times); end

fprintf('Default time is %f, if this is OK, <cr>;\n',time)
fprintf('If not, enter a time between %.4f and %.4f, we use the closest.\n', ...
                         min(times),max(times))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), time  = str2num(varstring); end 

d       = abs(time - times);    % crude distance
ind     = find(min(d) == d);  % multiple minima possible 
timeind = ind(1);             % use the first one
time    = times(timeind);


function [level, lvlind] = GetLevel(pgvar, levels, deflevel)
%----------------------------------------------------------------------
% level and lvlind will not be equal for all models, (and probably
% shouldn't for cam ... but for future expansion ...
if (nargin == 3), lvlind = deflevel; else lvlind = 1; end

if strcmpi(pgvar,'ps') ==1 
   disp('''PS'' only has one level, using it.')
   level  = 1;
   lvlind = 1;
else
   fprintf('Default level (index) is  %d, if this is OK, <cr>;\n',lvlind)
   fprintf('If not, enter a level between %d and %d, inclusive ...\n', ...
                         1,length(levels))
   varstring = input('we''ll use the closest (no syntax required)\n','s');

   if ~isempty(varstring), lvlind = str2num(varstring); end 

   % d      = abs(level - levels);  % crude distance
   % ind    = find(min(d) == d);    % multiple minima possible 
   % lvlind = ind(1);               % use the first one
   level  = levels(lvlind);
end



function [lon, lonind] = GetLongitude(pgvar, lons, deflon)
%----------------------------------------------------------------------
if (nargin == 3), lon = deflon; else lon = 255.0; end

fprintf('Default longitude is %f, if this is OK, <cr>;\n',lon)  
fprintf('If not, enter a longitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lons),max(lons))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lon  = str2num(varstring); end 

d      = abs(lon - lons);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible 
lonind = ind(1);             % use the first one
lon    = lons(lonind);



function [lat, latind] = GetLatitude(pgvar, lats, deflat)
%----------------------------------------------------------------------
if (nargin == 3), lat = deflat; else lat = 40.0; end

fprintf('Default latitude is %f, if this is OK, <cr>;\n',lat)
fprintf('If not, enter a latitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lats),max(lats))
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
