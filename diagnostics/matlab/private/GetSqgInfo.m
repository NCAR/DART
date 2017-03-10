function pinfo = GetSqgInfo(pstruct,fname,routine)
%% GetSqgInfo   prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetSqgInfo(pinfo_in,fname,routine);
%
% pinfo_in  Name of existing pinfo struct, e.g. output from CheckModelCompatibility
% fname     Name of the DART netcdf file
% routine   name of subsequent plot routine.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

pinfo  = pstruct;
model  = ncreadatt(fname,'/','model');

if strcmpi(model,'sqg') ~= 1
   error('Not so fast, this is not a sqg model.')
end

copy   = ncread(fname,'copy');
times  = ncread(fname,'time');
levels = ncread(fname,'lev');
lons   = ncread(fname,'lon');
lats   = ncread(fname,'lat');

switch lower(deblank(routine))

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar           = GetVar(pinfo.vars);
      [level, lvlind] = GetLevel(levels);
      [lat  , latind] = GetLatitude(lats);
      [lon  , lonind] = GetLongitude(lons);

      pinfo.var        = pgvar;
      pinfo.level      = level;
      pinfo.levelindex = lvlind;
      pinfo.longitude  = lon;
      pinfo.lonindex   = lonind;
      pinfo.latitude   = lat;
      pinfo.latindex   = latind;


   case 'plotcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(pinfo.vars);
      [base_time, base_tmeind] = GetTime(pinfo.time);
      [base_lvl,  base_lvlind] = GetLevel(levels);
      [base_lat,  base_latind] = GetLatitude(lats);
      [base_lon,  base_lonind] = GetLongitude(lons);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(pinfo.vars,         base_var);
      [comp_lvl, comp_lvlind] = GetLevel(levels, base_lvlind);

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
       base_var                = GetVar(pinfo.vars);
      [base_time, base_tmeind] = GetTime(pinfo.time);
      [base_lvl , base_lvlind] = GetLevel(levels);
      [base_lat , base_latind] = GetLatitude(lats);
      [base_lon , base_lonind] = GetLongitude(lons);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(pinfo.vars, base_var);
      [comp_lvl, comp_lvlind] = GetLevel(levels, base_lvlind);
      [comp_lat, comp_latind] = GetLatitude(lats,  base_lat);
      [comp_lon, comp_lonind] = GetLongitude(lons, base_lon);

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

       pgvar          = GetVar(pinfo.vars);
      [level, lvlind] = GetLevel(levels);
      [lat  , latind] = GetLatitude(lats);
      [lon  , lonind] = GetLongitude(lons);
      copyindices     = SetCopyID2(fname);
      copy            = length(copyindices);

      pinfo.var_names      = pgvar;
      pinfo.truth_file     = [];
      pinfo.prior_file     = pinfo.prior_file;
      pinfo.posterior_file = pinfo.posterior_file;
      pinfo.level          = level;
      pinfo.levelindex     = lvlind;
      pinfo.latitude       = lat;
      pinfo.latindex       = latind;
      pinfo.longitude      = lon;
      pinfo.lonindex       = lonind;
      pinfo.copies         = copy;
      pinfo.copyindices    = copyindices;

      if ( exist(pstruct.truth_file,'file') )
         pinfo.truth_file = pstruct.truth_file;
      end

   case 'plotphasespace'

      disp('Getting information for the ''X'' variable.')
       var1                   = GetVar(pinfo.vars);
      [var1_lvl, var1_lvlind] = GetLevel(levels);
      [var1_lat, var1_latind] = GetLatitude(lats);
      [var1_lon, var1_lonind] = GetLongitude(lons);

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVar(pinfo.vars, var1);
      [var2_lvl, var2_lvlind] = GetLevel(levels, var1_lvlind);
      [var2_lat, var2_latind] = GetLatitude(lats, var1_lat);
      [var2_lon, var2_lonind] = GetLongitude(lons, var1_lon);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVar(pinfo.vars, var1);
      [var3_lvl, var3_lvlind] = GetLevel(levels, var1_lvlind);
      [var3_lat, var3_latind] = GetLatitude(lats, var1_lat);
      [var3_lon, var3_lonind] = GetLongitude(lons, var1_lon);

      % query for copy string ... ensemble member
      [~, ens_mem_string] = GetCopyID(fname,1);
      ens_mem = ens_mem_string{1};   % coerce to simple character string

  %   s1 = input('Input ensemble member metadata ID. <cr> for ''true state''  ','s');
  %   if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''  ','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

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

   case 'plottotalerr'

     % Nothing required ... function defined in PlotTotalErr()

   otherwise

      error('plot type %s not supported yet for SQG model',routine)

end



function pgvar = GetVar(prognostic_vars, defvar)
%----------------------------------------------------------------------
if (nargin == 2), pgvar = defvar; else pgvar = 'theta'; end

% If there is only one choice ... use it.
if (length(prognostic_vars) == 1)
   pgvar = prognostic_vars{1};
   return
end

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



function [time, timeind] = GetTime(times, deftime)
%----------------------------------------------------------------------
% Query for the time of interest.

ntimes = length(times);

% Determine a sensible default.
if (nargin == 2),
   time = deftime;
   tindex = find(times == deftime);
else
   if (ntimes < 2)
      tindex = round(ntimes/2);
   else
      tindex = 1;
   end
   time = times(tindex);
end

fprintf('Default time is %s (index %d), if this is OK, <cr>;\n',datestr(time),tindex)
fprintf('If not, enter an index between %d and %d \n',1,ntimes)
fprintf('Pertaining to %s and %s \n',datestr(times(1)),datestr(times(ntimes)))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), tindex = str2double(varstring); end

timeinds = 1:ntimes;
d        = abs(tindex - timeinds); % crude distance
ind      = find(min(d) == d);      % multiple minima possible
timeind  = ind(1);                 % use the first one
time     = times(timeind);


function [level, lvlind] = GetLevel(levels, deflevel)
%----------------------------------------------------------------------
% level and lvlind will not be equal ... future expansion ...
if (nargin == 3), lvlind = deflevel; else lvlind = 1; end

fprintf('Default level (index) is  %d, if this is OK, <cr>;\n',lvlind)
fprintf('If not, enter a level between %d and %d, inclusive ...\n', ...
                     1,length(levels))
varstring = input('we''ll use the closest (no syntax required)\n','s');

if ~isempty(varstring), lvlind = str2double(varstring); end

level  = levels(lvlind);



function [lon, lonind] = GetLongitude(lons, deflon)
%----------------------------------------------------------------------
if (nargin == 2), lon = deflon; else lon = 255.0; end

fprintf('Default longitude is %f, if this is OK, <cr>;\n',lon)
fprintf('If not, enter a longitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lons),max(lons))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lon  = str2double(varstring); end

d      = abs(lon - lons);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
lonind = ind(1);             % use the first one
lon    = lons(lonind);



function [lat, latind] = GetLatitude(lats, deflat)
%----------------------------------------------------------------------
if (nargin == 2), lat = deflat; else lat = 40.0; end

fprintf('Default latitude is %f, if this is OK, <cr>;\n',lat)
fprintf('If not, enter a latitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lats),max(lats))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lat = str2double(varstring); end

d      = abs(lat - lats);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
latind = ind(1);             % use the first one
lat    = lats(latind);




function [copyid, copystrings] = GetCopyID(fname, numcopies)
%----------------------------------------------------------------------
%% GetCopyID queries for the copy indices in a specific netCDF file.

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

metadata = nc_varget(fname,'CopyMetaData');           % get all the metadata
ncopies  = size(metadata,1);

if ( ncopies < 1 )
   fprintf('%s has no valid ensemble members\n',fname)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',fname, numcopies)
   error('netcdf file has no ensemble members.')
end

if (numcopies > ncopies)
   error('Only have %d copies to work with, need %d',ncopies,numcopies)
end

for i = 1:ncopies
   fprintf('ID %2d  is  %s\n',i,deblank(metadata(i,:)))
end

IDstring = input( sprintf('Enter %d ID(s) to plot.\n(no intervening syntax, please)\n',numcopies) ,'s');
bob      = str2num(IDstring);

% If entered more than desired, just use first N
if (length(bob) >= numcopies)
   copyid = bob(1:numcopies);
else
   error('only entered %d, needed %d ... quitting',length(bob),numcopies)
end

copystrings = cell(1,numcopies);
for i=1:numcopies
  copystrings{i} = metadata(copyid(i),:);
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
