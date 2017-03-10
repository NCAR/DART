function pinfo = GetPOPInfo(pinfo_in,fname,routine)
%% GetPOPInfo   prepares a structure of information needed by the subsequent "routine"
%               The information is gathered via rudimentary "input" routines.
%
% pinfo = GetPOPInfo(pinfo_in, fname, routine);
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

pinfo = pinfo_in;

if strcmpi(pinfo.model,'POP') ~= 1
    error('Not so fast, this is not a POP model.')
end

%% Get the required information.

global copy ULON ULAT TLON TLAT ZG ZC KMU KMT

varexist(fname, {'copy','ULON','ULAT','TLON','TLAT','ZG','ZC','KMT','KMU'});

copy = ncread(fname, 'copy');
ULON = ncread(fname, 'ULON');
ULAT = ncread(fname, 'ULAT');
TLON = ncread(fname, 'TLON');
TLAT = ncread(fname, 'TLAT');
ZG   = ncread(fname,   'ZG');
ZC   = ncread(fname,   'ZC');
KMU  = ncread(fname,  'KMU');
KMT  = ncread(fname,  'KMT');

depthunits = ncreadatt(fname,'ZG','units');


%% The POP model has different number of possible levels for each location.
%  Consequently, it is required to query the lat/lon before the number of
%  levels.

switch lower(deblank(routine))

    case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

        pgvar           = GetVarString(pinfo_in.vars);     % set prognostic variable
        [lat, lon, latind, lonind] = GetLatLon(fname, pgvar);
        [level, lvlind] = GetLevel(fname,pgvar,lonind,latind);

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
        base_var                 = GetVarString(pinfo_in.vars);
        [base_time, base_tmeind] = GetTime(pinfo.time);
        [base_lat, base_lon, base_latind, base_lonind] = GetLatLon(fname, base_var);
        [base_lvl,  base_lvlind] = GetLevel(fname, base_var, base_lonind, base_latind);

        disp('Getting information for the ''comparison'' variable.')
        comp_var                = GetVarString(pinfo_in.vars, base_var);
        [comp_lvl, comp_lvlind] = GetLevel(fname, base_var, base_lonind, base_latind, base_lvl);

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
        pinfo.depthunits  = depthunits;

    case 'plotvarvarcorrel'

        disp('Getting information for the ''base'' variable.')
        base_var                 = GetVarString(pinfo_in.vars);
        [base_time, base_tmeind] = GetTime(pinfo.time);
        [base_lat, base_lon, base_latind, base_lonind] = GetLatLon(fname, base_var);
        [base_lvl , base_lvlind] = GetLevel(fname, base_var, base_lonind, base_latind);

        disp('Getting information for the ''comparison'' variable.')
        comp_var                = GetVarString(pinfo_in.vars, base_var);
        [comp_lat, comp_lon, comp_latind, comp_lonind] = GetLatLon(fname, comp_var, ...
                                                            base_latind, base_lonind);
        [comp_lvl, comp_lvlind] = GetLevel(fname, comp_var, comp_lonind, comp_latind, base_lvl);

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
        pinfo.depthunits  = depthunits;

    case 'plotsawtooth'

        pgvar                      = GetVarString(pinfo_in.vars);
        [lat, lon, latind, lonind] = GetLatLon(pinfo_in.prior_file, pgvar);
        [level, lvlind]            = GetLevel(pinfo_in.prior_file, pgvar, lonind, latind);
        [copyindices, copymetadata]= SetCopyID2(pinfo_in.prior_file);
        copy                       = length(copyindices);

        pinfo.truth_file     = pinfo_in.truth_file;
        pinfo.prior_file     = pinfo_in.prior_file;
        pinfo.posterior_file = pinfo_in.posterior_file;
        pinfo.var_names      = pgvar;
        pinfo.level          = level;
        pinfo.levelindex     = lvlind;
        pinfo.latitude       = lat;
        pinfo.latindex       = latind;
        pinfo.longitude      = lon;
        pinfo.lonindex       = lonind;
        pinfo.copies         = copy;
        pinfo.copyindices    = copyindices;
        pinfo.copymetadata   = copymetadata;

        % XXX has a habit of cutting out the rest of the ensemble members
        % from the posterior, but leaving them in the prior. Thanks.
        % So now I have to figure out if the posterior and prior copy metadata match.

        for i = 1:copy,
            copyi = get_member_index(pinfo_in.posterior_file,copymetadata{i});
            pstruct.postcopyindices = copyi;
        end

        if ( any(pstruct.postcopyindices < 1) )
            error('The prior copy does not exist in the posterior file.')
        end

    case 'plotphasespace'

        disp('Getting information for the ''X'' variable.')
        var1                   = GetVarString(pinfo_in.vars);
        [var1_lat, var1_lon, var1_latind, var1_lonind] = GetLatLon(fname, var1);
        [var1_lvl, var1_lvlind] = GetLevel(fname,var1,var1_lonind,var1_latind);

        disp('Getting information for the ''Y'' variable.')
        var2                   = GetVarString(pinfo_in.vars,        var1    );
        [var2_lat, var2_lon, var2_latind, var2_lonind] = GetLatLon(fname, var2, var1_latind, var1_lonind);
        [var2_lvl, var2_lvlind] = GetLevel(fname,var2,var2_lonind,var2_latind,var1_lvl);

        disp('Getting information for the ''Z'' variable.')
        var3                   = GetVarString(pinfo_in.vars,        var1    );
        [var3_lat, var3_lon, var3_latind, var3_lonind] = GetLatLon(fname, var3, var1_latind, var1_lonind);
        [var3_lvl, var3_lvlind] = GetLevel(fname,var3,var3_lonind,var3_latind,var1_lvl);

        % query for ensemble member string
        metadata   = ncread(fname,'MemberMetadata');
        [N,M]      = size(metadata);
        if M == 1
            cell_array{1} = metadata';
        else
            cell_array = mat2cell(metadata, ones(1,N), M);
        end
        ens_mem    = strtrim(cell_array{1});
        str1 = sprintf('Input ensemble member metadata STRING. <cr> for ''%s''   ',ens_mem);
        s1   = input(str1,'s');
        if ~ isempty(s1), ens_mem = s1; end

        % query for line type
        s1 = input('Input line type string. <cr> for ''k-''  ','s');
        if isempty(s1), ltype = 'k-'; else ltype = s1; end

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

if ~isempty(varstring), tindex = str2num(varstring); end

timeinds = 1:ntimes;
d        = abs(tindex - timeinds); % crude distance
ind      = find(min(d) == d);      % multiple minima possible
timeind  = ind(1);                 % use the first one
time     = times(timeind);



function [level, lvlind] = GetLevel(fname, pgvar, lonind, latind, deflevel)
%----------------------------------------------------------------------
%
global ZC ZG KMU KMT

if (nargin == 5), level = deflevel; else level = []; end

varinfo  = nc_getvarinfo(fname,pgvar);
leveldim = find(strncmp('k',varinfo.Dimension,1) > 0);

if (isempty(leveldim))

    fprintf('''%s'' has no vertical level, imply using 1.\n',pgvar)
    level  = 1;
    lvlind = 1;

else

    levelvar    = varinfo.Dimension{leveldim};
    [nlevels,~] = nc_dim_info(fname,levelvar);
    levels      = 1:nlevels;

    % must find the lowest level at this gridcell for this variable
    velcomp = regexp(varinfo.Name,'VEL','ONCE');

    if ( isempty(velcomp) )
        levbottom = KMT(lonind,latind);
        levels    = ZC(1:levbottom);
    else
        levbottom = KMU(lonind,latind);
        levels    = ZG(1:levbottom);
    end

    if (isempty(level)), level = levels(1); end

    fprintf('Default level is  %0.2f, if this is OK, <cr>;\n',level)
    fprintf('If not, enter a level between %0.2f and %0.2f, inclusive ...\n', ...
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
% The POP domain is not a simple 1D lat 1D lon domain, so it makes more
% sense to query about lat & lon indices rather than values.

global ULON ULAT TLON TLAT

varinfo = nc_getvarinfo(fname, pgvar);

londim  = find(strncmp('i',varinfo.Dimension,1) > 0);
latdim  = find(strncmp('j',varinfo.Dimension,1) > 0);

if (isempty(latdim) || isempty(londim))
    error('''%s'' does not have both lats and lons.\n',pgvar)
end

nlat = varinfo.Size(latdim);
nlon = varinfo.Size(londim);

% both U,V velocity components are on the ULON, ULAT coordinates.
% The T,S components use TLON, TLAT

velcomp = regexp(varinfo.Name,'VEL','ONCE');

if ( isempty(velcomp) )
    latmat = TLAT;
    lonmat = TLON;
else
    latmat = ULAT;
    lonmat = ULON;
end

% Determine a sensible default.
if (nargin > 2)
    latind = deflat;
    lonind = deflon;
else % in the middle of the domain.
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
        error('Did not get two indices for the lat lon pair.')
    end
    latind = nums(1);
    lonind = nums(2);
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

d      = abs(latind - latinds);  % crude distance
ind    = find(min(d) == d);   % multiple minima possible
latind = ind(1);              % use the first one

d      = abs(lonind - loninds);  % crude distance
ind    = find(min(d) == d);   % multiple minima possible
lonind = ind(1);              % use the first one

lat    = latmat(latind,lonind);
lon    = lonmat(latind,lonind);



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



% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
