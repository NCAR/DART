function pinfo = GetMITgcm_oceanInfo(pstruct,fname,routine);
%% GetMITgcm_oceanInfo   prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetMITgcm_oceanInfo(fname,routine);
%
% fname     Name of the DART netcdf file
% routine   name of subsequent plot routine.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

if strcmp(lower(pstruct.model),'mitgcm_ocean') ~= 1
   error('Not so fast, this is not an MITgcm_ocean model.')
end

copy   = ncread(fname,'copy');
levels = ncread(fname,  'ZG');
XG     = ncread(fname,  'XG');
XC     = ncread(fname,  'XC');
YG     = ncread(fname,  'YG');
YC     = ncread(fname,  'YC');

times      = ncread(fname,'time');
timeunits  = ncreadatt(fname,'time','units');
timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
dates      = times + timeorigin;
num_times  = length(dates);

disp('GetMITgcm_oceanInfo: pstruct should have time/dates already ...')

clear times timeunits timebase timeorigin


% A more robust way would be to use the netcdf low-level ops:
% bob = var(f);     % bob is a cell array of ncvars
% name(bob{1})       % is the variable name string
% bob{1}(:)          % is the value of the netcdf variable  (no
% offset/scale)
% have not yet figured out a way to only use non-coordinate variables.

prognostic_vars = {'S','T','U','V','SSH'};
num_vars = length(prognostic_vars);

pinfo = pstruct;

switch lower(deblank(routine))

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar           = GetVar(prognostic_vars);  % Determine prognostic variable
      [level, lvlind] = GetLevel(pgvar,levels);   % Determine level and index
      [lat  , latind] = GetLatitude( pgvar,YG,YC);
      [lon  , lonind] = GetLongitude(pgvar,XG,XC);

      pinfo.fname      = fname;
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
      [base_lat,  base_latind] = GetLatitude( base_var,YG,YC);
      [base_lon,  base_lonind] = GetLongitude(base_var,XG,XC);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvl);

      pinfo.fname     = fname      ;
      pinfo.base_var  = base_var ; pinfo.comp_var    = comp_var   ;
      pinfo.base_time = base_time; pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl  = base_lvl ; pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat  = base_lat ; pinfo.base_latind = base_latind;
      pinfo.base_lon  = base_lon ; pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl  = comp_lvl ; pinfo.comp_lvlind = comp_lvlind;

   case 'plotvarvarcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(prognostic_vars);
      [base_time, base_tmeind] = GetTime(     base_var,times);
      [base_lvl , base_lvlind] = GetLevel(    base_var,levels);
      [base_lat , base_latind] = GetLatitude( base_var,YG,YC);
      [base_lon , base_lonind] = GetLongitude(base_var,XG,XC);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvl);
      [comp_lat, comp_latind] = GetLatitude( comp_var,YG,YC);
      [comp_lon, comp_lonind] = GetLongitude(comp_var,XG,XC);

      pinfo.fname     = fname      ;
      pinfo.base_var  = base_var ; pinfo.comp_var    = comp_var   ;
      pinfo.base_time = base_time; pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl  = base_lvl ; pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat  = base_lat ; pinfo.base_latind = base_latind;
      pinfo.base_lon  = base_lon ; pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl  = comp_lvl ; pinfo.comp_lvlind = comp_lvlind;
      pinfo.comp_lat  = comp_lat ; pinfo.comp_latind = comp_latind;
      pinfo.comp_lon  = comp_lon ; pinfo.comp_lonind = comp_lonind;

   case 'plotsawtooth'

       pgvar          = GetVar(prognostic_vars);
      [level, lvlind] = GetLevel(    pgvar,levels);
      [  lat, latind] = GetLatitude( pgvar,YG,YC);
      [  lon, lonind] = GetLongitude(pgvar,XG,XC);
      copyindices     = SetCopyID(fname);
      copy            = length(copyindices);

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

   case 'plotphasespace'

      disp('Getting information for the ''X'' variable.')
       var1                   = GetVar(prognostic_vars);
      [var1_lvl, var1_lvlind] = GetLevel(    var1, levels);
      [var1_lat, var1_latind] = GetLatitude( var1, YG, YC);
      [var1_lon, var1_lonind] = GetLongitude(var1, XG, XC );

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVar(prognostic_vars,    var1    );
      [var2_lvl, var2_lvlind] = GetLevel(    var2, levels, var1_lvl);
      [var2_lat, var2_latind] = GetLatitude( var2, YG, YC, var1_lat);
      [var2_lon, var2_lonind] = GetLongitude(var2, XG, XC, var1_lon);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVar(prognostic_vars,    var1    );
      [var3_lvl, var3_lvlind] = GetLevel(    var3, levels, var1_lvl);
      [var3_lat, var3_latind] = GetLatitude( var3, YG, YC, var1_lat);
      [var3_lon, var3_lonind] = GetLongitude(var3, XG, XC, var1_lon);

      % query for ensemble member
      s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
      if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''  ','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

      pinfo.fname    =fname;
      pinfo.var1name = var1;
      pinfo.var2name = var2;
      pinfo.var3name = var3;
      pinfo.var1_lvl = var1_lvl; pinfo.var1_lvlind = var1_lvlind;
      pinfo.var1_lat = var1_lat; pinfo.var1_latind = var1_latind;
      pinfo.var1_lon = var1_lon; pinfo.var1_lonind = var1_lonind;
      pinfo.var2_lvl = var2_lvl; pinfo.var2_lvlind = var2_lvlind;
      pinfo.var2_lat = var2_lat; pinfo.var2_latind = var2_latind;
      pinfo.var2_lon = var2_lon; pinfo.var2_lonind = var2_lonind;
      pinfo.var3_lvl = var3_lvl; pinfo.var3_lvlind = var3_lvlind;
      pinfo.var3_lat = var3_lat; pinfo.var3_latind = var3_latind;
      pinfo.var3_lon = var3_lon; pinfo.var3_lonind = var3_lonind;
      pinfo.ens_mem  = ens_mem ; pinfo.ltype       = ltype;

   otherwise

end

function pgvar = GetVar(prognostic_vars, defvar)
%----------------------------------------------------------------------
if (nargin == 2), pgvar = defvar; else pgvar = 'S'; end

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
a           = strread(pgvar,'%s','delimiter',' ');


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
% shouldn't for the bgrid ... but for future expansion ...
if (nargin == 3)
    lvlind = find(levels == deflevel);
else
    lvlind = 1;
end

if strcmp(lower(pgvar),'ssh') ==1
   disp('''SSH'' only has one level, using it.')
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



function [lon, lonind] = GetLongitude(pgvar, xg, xc, deflon)
%----------------------------------------------------------------------
if (nargin == 4), lon = deflon; else lon = 270.0; end

switch lower(pgvar)
   case {'u'}
     lons = xg;
   otherwise
     lons = xc;
end

fprintf('Default longitude is %f, if this is OK, <cr>;\n',lon)
fprintf('If not, enter a longitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lons),max(lons))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lon  = str2num(varstring); end

d      = abs(lon - lons);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
lonind = ind(1);             % use the first one
lon    = lons(lonind);



function [lat, latind] = GetLatitude(pgvar, yg, yc, deflat)
%----------------------------------------------------------------------
if (nargin == 4), lat = deflat; else lat = 25.0; end

switch lower(pgvar)
   case {'v'}
     lats = yg;
   otherwise
     lats = yc;
end

fprintf('Default latitude is %f, if this is OK, <cr>;\n',lat)
fprintf('If not, enter a latitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lats),max(lats))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lat = str2num(varstring); end

d      = abs(lat - lats);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
latind = ind(1);             % use the first one
lat    = lats(latind);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
