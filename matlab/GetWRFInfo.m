function pinfo = GetWRFInfo(pinfo_in,fname,routine)
%% GetWRFInfo   prepares a structure of information needed by the subsequent "routine"
%               The information is gathered via rudimentary "input" routines.
%
% pinfo = GetWRFInfo(pinfo_in, fname, routine);
%
% pinfo_in  Name of existing pinfo struct, e.g. output from CheckModelCompatibility
% fname     Name of the DART netcdf file
% routine   name of subsequent plot routine.
%
% fname = '/glade/proj2/image/romine/dart/work_Radar/rad_regression/geom/Prior_Diag.nc';
%

%% DART software - Copyright � 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

pinfo = pinfo_in;
model = nc_attget(fname, nc_global, 'model');

if strcmpi(model,'wrf') ~= 1
   error('Not so fast, this is not a WRF model.')
end

%% Get the domain-independent information.

varexist(fname, {'copy','time'})

copy       = nc_varget(fname,'copy');
times      = nc_varget(fname,'time');

dx         = varget(fname,        'DX');
dy         = varget(fname,        'DY');
truelat1   = varget(fname,  'TRUELAT1');
truelat2   = varget(fname,  'TRUELAT2');
stand_lon  = varget(fname, 'STAND_LON');
cen_lat    = varget(fname,   'CEN_LAT');
cen_lon    = varget(fname,   'CEN_LON');
map_proj   = varget(fname,  'MAP_PROJ');
periodic_x = varget(fname,'PERIODIC_X');
polar      = varget(fname,     'POLAR');

%% Get the global metadata from a WRF DART diagnostic netCDF file.
% If there is only one domain, we know what to do. 
% otherwise, ask which domain is of interest.

dinfo       = nc_getdiminfo(fname,'domain');  % no graceful error
num_domains = dinfo.Length;

dID    = 1;
if (num_domains > 1) 
   dID = GetDomain(num_domains);  % Determine prognostic variable
end

dn      = varget(fname, sprintf(     'DN_d%02d',dID));
znu     = varget(fname, sprintf(    'ZNU_d%02d',dID));
znw     = varget(fname, sprintf(    'ZNW_d%02d',dID));
dnw     = varget(fname, sprintf(    'DNW_d%02d',dID));
mub     = varget(fname, sprintf(    'MUB_d%02d',dID));
xlong   = varget(fname, sprintf(  'XLONG_d%02d',dID));
xlong_u = varget(fname, sprintf('XLONG_U_d%02d',dID));
xlong_v = varget(fname, sprintf('XLONG_V_d%02d',dID));
xlat    = varget(fname, sprintf(   'XLAT_d%02d',dID));
xlat_u  = varget(fname, sprintf( 'XLAT_U_d%02d',dID));
xlat_v  = varget(fname, sprintf( 'XLAT_V_d%02d',dID));
levels  = varget(fname, sprintf(  'level_d%02d',dID));
xland   = varget(fname, sprintf(  'XLAND_d%02d',dID));
phb     = varget(fname, sprintf(    'PHB_d%02d',dID));
hgt     = varget(fname, sprintf(    'HGT_d%02d',dID));

prognostic_bases = {'U', 'V', 'W', 'PH', 'T', ...
          'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE'};

for i = 1:length(prognostic_bases)
   prognostic_vars{i} = sprintf('%s_d%02d',prognostic_bases{i},dID);
end

varexist(fname, prognostic_vars)

switch lower(deblank(routine))

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar           = GetVarString(prognostic_vars); % Determine prognostic variable
      [level, lvlind] = GetLevel(fname,pgvar);         % Determine level and index
      [lat, lon, latind, lonind] = GetLatLon(fname, pgvar);

      pinfo.model      = model;
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
       base_var                = GetVarString(prognostic_vars);
      [base_time, base_tmeind] = GetTime(times);
      [base_lvl,  base_lvlind] = GetLevel(fname,base_var);
      [base_lat,  base_latind] = GetLatitude( base_var,TmpJ,VelJ);
      [base_lon,  base_lonind] = GetLongitude(base_var,TmpI,VelI);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVarString(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(fname,comp_var,base_lvl);

      pinfo.model       = model;
      pinfo.fname       = fname;
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
       base_var                = GetVarString(prognostic_vars);
      [base_time, base_tmeind] = GetTime(times);
      [base_lvl , base_lvlind] = GetLevel(fname,base_var);
      [base_lat , base_latind] = GetLatitude( base_var,TmpJ,VelJ);
      [base_lon , base_lonind] = GetLongitude(base_var,TmpI,VelI);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVarString(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(fname,comp_var,base_lvl);
      [comp_lat, comp_latind] = GetLatitude( comp_var,TmpJ,VelJ, base_lat);
      [comp_lon, comp_lonind] = GetLongitude(comp_var,TmpI,VelI, base_lon);

      pinfo.model       = model;
      pinfo.fname       = fname;
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

       pgvar          = GetVarString(prognostic_vars);
      [level, lvlind] = GetLevel(fname,pgvar);   % Determine level and index
      [lat  , latind] = GetLatitude( pgvar,TmpJ,VelJ);
      [lon  , lonind] = GetLongitude(pgvar,TmpI,VelI);
      %[copy , lonind] = GetCopies(pgvar,copy);
      copyindices     = SetCopyID(fname);
      copy            = length(copyindices);

      pinfo.model          = model;
      pinfo.var_names      = pgvar;
     %pinfo.truth_file     = [];
     %pinfo.prior_file     = pinfo.prior_file;
     %pinfo.posterior_file = pinfo.posterior_file;
      pinfo.level          = level;
      pinfo.levelindex     = lvlind;
      pinfo.latitude       = lat;
      pinfo.latindex       = latind;
      pinfo.longitude      = lon;
      pinfo.lonindex       = lonind;
      pinfo.copies         = copy;
      pinfo.copyindices    = copyindices;

   case 'plotphasespace'

      disp('Getting information for the ''X'' variable.')
       var1                   = GetVarString(prognostic_vars);
      [var1_lvl, var1_lvlind] = GetLevel(fname,var1);
      [var1_lat, var1_latind] = GetLatitude( var1, TmpJ, VelJ);
      [var1_lon, var1_lonind] = GetLongitude(var1, TmpI, VelI);

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVarString(prognostic_vars,        var1    );
      [var2_lvl, var2_lvlind] = GetLevel(fname,var2,var1_lvl);
      [var2_lat, var2_latind] = GetLatitude( var2, TmpJ, VelJ, var1_lat);
      [var2_lon, var2_lonind] = GetLongitude(var2, TmpI, VelI, var1_lon);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVarString(prognostic_vars,        var1    );
      [var3_lvl, var3_lvlind] = GetLevel(fname,var3,var1_lvl);
      [var3_lat, var3_latind] = GetLatitude( var3, TmpJ, VelJ, var1_lat);
      [var3_lon, var3_lonind] = GetLongitude(var3, TmpI, VelI, var1_lon);

      % query for ensemble member
      s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
      if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''  ','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

      pinfo.model       = model;
      pinfo.fname       = fname;
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

function domainid = GetDomainID(ndomains)
%----------------------------------------------------------------------

fprintf('There are %d domains. Default is to use domain %d.\n.',ndomains,ndomains)  
fprintf('If this is OK, <cr>; If not, enter domain of interest:\n')
domainid = input('(no syntax required)\n');

if ~isempty(domainid), domainid = ndomains; end 

if ( (domainid > 0) && (domainid <= ndomains)) 
   error('domain must be between 1 and %d, you entered %d',ndomains,domainid)
end


function pgvar = GetVarString(prognostic_vars, defvar)
%----------------------------------------------------------------------
pgvar = prognostic_vars{1};
if (nargin == 2), pgvar = defvar; end

str = sprintf(' %s ',prognostic_vars{1});
for i = 2:length(prognostic_vars),
   str = sprintf(' %s %s ',str,prognostic_vars{i});
end
fprintf('Default variable is ''%s'', if this is OK, <cr>;\n',pgvar)  
fprintf('If not, please enter one of: %s\n',str)
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), pgvar = deblank(varstring); end 



function [time, timeind] = GetTime(times, deftime)
%----------------------------------------------------------------------
if (nargin == 2), time = deftime; else time = mean(times); end

fprintf('Default time is %f, if this is OK, <cr>;\n',time)
fprintf('If not, enter a time between %.4f and %.4f, we use the closest.\n', ...
                         min(times),max(times))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), time  = str2num(varstring); end 

d       = abs(time - times);  % crude distance
ind     = find(min(d) == d);  % multiple minima possible 
timeind = ind(1);             % use the first one
time    = times(timeind);


function [level, lvlind] = GetLevel(fname, pgvar, deflevel)
%----------------------------------------------------------------------
%
if (nargin == 3), level = deflevel; else level = []; end

varinfo  = nc_getvarinfo(fname,pgvar);
leveldim = find(strncmp('bottom_top',varinfo.Dimension,length('bottom_top')) > 0);

if (isempty(leveldim))

   fprintf('''%s'' has no vertical level, imply using 1.\n',pgvar)
   level  = 1;
   lvlind = 1;

else

   levelvar = varinfo.Dimension{leveldim};
   dinfo    = nc_getdiminfo(fname,levelvar); 
   levels   = 1:dinfo.Length;

   if (isempty(level)), level = levels(1); end

   fprintf('Default level is  %d, if this is OK, <cr>;\n',level)
   fprintf('If not, enter a level between %d and %d, inclusive ...\n', ...
                         min(levels),max(levels))
   varstring = input('we''ll use the closest (no syntax required)\n','s');

   if ~isempty(varstring), level = str2num(varstring); end 

   d      = abs(level - levels);  % crude distance
   ind    = find(min(d) == d);    % multiple minima possible 
   lvlind = ind(1);               % use the first one
   level  = levels(lvlind);

end



function [lat, lon, latind, lonind] = GetLatLon(fname, pgvar, deflat, deflon)
%----------------------------------------------------------------------
% The WRF domain is not a simple 1D lat 1D lon domain, so it makes more
% sense to query about lat & lon indices rather than values.

varinfo = nc_getvarinfo(fname, pgvar);
latdim  = find(strncmp('south_north',varinfo.Dimension,length('south_north')) > 0);
londim  = find(strncmp(  'west_east',varinfo.Dimension,length(  'west_east')) > 0);
nlat    = varinfo.Size(latdim);
nlon    = varinfo.Size(londim);

if (isempty(latdim) || isempty(londim))
   error('''%s'' does not have both lats and lons.\n',pgvar)
end

% Are either of the grids staggered. isempty() means unstaggered.
latstag = regexp(varinfo.Dimension{latdim},'.stag.');
lonstag = regexp(varinfo.Dimension{londim},'.stag.');

% get the domain id from the pgvar string

dID = 1;

if ( isempty(latstag) && isempty(lonstag) )

   % Both unstaggered.
   latmat = nc_varget(fname, 'XLAT_d01');
   lonmat = nc_varget(fname,'XLONG_d01');

elseif ( isempty(latstag) )

   % LAT unstaggered, LON staggered.
   latmat = nc_varget(fname, 'XLAT_U_d01');
   lonmat = nc_varget(fname,'XLONG_U_d01');

else

   % LAT staggered, LON unstaggered.
   latmat = nc_varget(fname, 'XLAT_V_d01');
   lonmat = nc_varget(fname,'XLONG_V_d01');

end


% Determine a sensible default.
if (nargin > 2)
   latind = deflat;
   lonind = deflon;
else
   latind = round(nlat/2);
   lonind = round(nlon/2);
end
lat = latmat(latind,lonind);
lon = lonmat(latind,lonind);

% Ask if they like the default.

fprintf('There are %d latitudes and %d longitudes for %s\n',nlat,nlon,pgvar)
fprintf('Default lat/lon indices are %d,%d ... lat/lon(%f,%f)\n',latind,lonind,lat,lon)
fprintf('If this is OK, <cr>; if not, enter grid indices lat/lon pair ...\n')
varstring = input('we''ll use the closest (no syntax required)\n','s');

if ~isempty(varstring)
   nums = str2num(varstring);
   if (length(nums) ~= 2) 
      error('Did not get two numbers for the lat lon pair.')
   end
   lat = nums(1);
   lon = nums(2);
end

% Find the closest lat/lon location to the requested one.
% not-so-Simple search over a matrix of lats/lons.

% londiff          = lonmat - lon;
% latdiff          = latmat - lat;
% dist             = sqrt(londiff.^2 + latdiff.^2);
% [latind, lonind] = find(min(dist(:)) == dist);
% lat              = latmat(latind,lonind);
% lon              = lonmat(latind,lonind);

% Much simpler search to find closest lat index.
latinds = 1:nlat;
loninds = 1:nlon;

d      = abs(lat - latinds);  % crude distance
ind    = find(min(d) == d);   % multiple minima possible 
latind = ind(1);              % use the first one

d      = abs(lon - loninds);  % crude distance
ind    = find(min(d) == d);   % multiple minima possible 
lonind = ind(1);              % use the first one

lat    = latmat(latind,lonind);
lon    = lonmat(latind,lonind);


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


function varexist(filename, varnames)
%% We already know the file exists by this point.
% Lets check to make sure that file contains all needed variables.
% Fatally die if they do not exist.

nvars  = length(varnames);
gotone = ones(1,nvars);

for i = 1:nvars
   gotone(i) = nc_isvar(filename,varnames{i});
   if ( ~ gotone(i) )
      fprintf('\n%s is not a variable in %s\n',varnames{i},filename)
   end
end

if ~ all(gotone)
   error('missing required variable ... exiting')
end



function x = varget(filename,varname)
%% get a varible from the netcdf file, if it does not exist, 
% do not die such a theatrical death ... return an empty.

if ( nc_isvar(filename,varname) )
   x = nc_varget(filename, varname);
else
   x = [];
end
