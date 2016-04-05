function pinfo = GetGCOM_oceanInfo(pstruct, fname, routine)
%% GetGCOM_oceanInfo  prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetGCOM_oceanInfo(fname,routine);
%
% fname     Name of the DART netcdf file
% routine   name of subsequent plot routine.

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 )
    error('%s does not exist.',fname)
end

if strcmpi(pstruct.model,'gcom') ~= 1
    error('Not so fast, this is not an GCOM_ocean model.')
end

% These will be different for each variable in GCOM
levels = nc_varget(fname, 'ulev');
XC     = nc_varget(fname, 'lon');
XG     = nc_varget(fname, 'ulon');
YC     = nc_varget(fname, 'lat');
YG     = nc_varget(fname, 'ulat');

% A more robust way would be to use the netcdf low-level ops:
% bob = var(f);     % bob is a cell array of ncvars
% name(bob{1})       % is the variable name string
% bob{1}(:)          % is the value of the netcdf variable  (no
% offset/scale)
% have not yet figured out a way to only use non-coordinate variables.

pinfo = pstruct;

% all the functions need at least one variable.
disp('Getting information for the ''base'' variable.')
metadata = getVariableProperties();

switch lower(deblank(routine))
    
    case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}
        
        % The fully curvilinear grid means the best way to query may be by
        % index and decode the lon/lat/lev at the end.
        
        lvlind = getLevelIndex();
        latind = getLatitudeIndex();
        lonind = getLongitudeIndex();
        
        pinfo.fname      = fname;
        pinfo.var        = metadata.varname;
        pinfo.lonindex   = lonind;
        pinfo.latindex   = latind;
        pinfo.levelindex = lvlind;
        pinfo.longitude  = metadata.longitude(lonind, latind, lvlind);
        pinfo.latitude   = metadata.latitude( lonind, latind, lvlind);
        pinfo.level      = metadata.level(    lonind, latind, lvlind);
        
    case 'plotcorrel'
        
        disp('Getting information for the ''base'' variable.')
        [base_time, base_tmeind] = GetTime(     base_var,times);
        base_lvlind = getLevelIndex();
        base_latind = getLatitudeIndex();
        base_lonind = getLongitudeIndex();
        
        disp('Getting information for the ''comparison'' variable.')
        comp_var_metadata = getVariableProperties(metadata);
        comp_lvlind       = getLevelIndex(base_lvlind);
        
        pinfo.fname       = fname;
        pinfo.base_var    = metadata.varname;
        pinfo.comp_var    = comp_var_metadata.varname;
        pinfo.base_tmeind = base_tmeind;
        pinfo.base_time   = metadata.time(base_tmeind);
        pinfo.base_lonind = base_lonind;
        pinfo.base_latind = base_latind;
        pinfo.base_lvlind = base_lvlind;
        pinfo.base_lon    = metadata.longitude(base_lonind, base_latind, base_lvlind);
        pinfo.base_lat    = metadata.latitude( base_lonind, base_latind, base_lvlind);
        pinfo.base_lvl    = metadata.level(    base_lonind, base_latind, base_lvlind);

        % The level is just for annotation purposes for the entire 'level' - 
        % cannot be correct for a general curvilinear coordinate system, but
        % we need to annotate something.
        pinfo.comp_lvlind = comp_lvlind;
        pinfo.comp_lvl    = comp_var_metadata.level(base_lonind, base_latind, comp_lvlind);
        
    case 'plotvarvarcorrel'
        
        disp('Getting information for the ''base'' variable.')
        [base_time, base_tmeind] = GetTime(     base_var,times,deftime);
        [base_lvl , base_lvlind] = getLevelIndex(    base_var,levels);
        [base_lat , base_latind] = getLatitudeIndex( base_var,YG,YC);
        [base_lon , base_lonind] = getLongitudeIndex(base_var,XG,XC);
        
        disp('Getting information for the ''comparison'' variable.')
        comp_var                = getVariableProperties(pinfo.vars,          base_var);
        [comp_lvl, comp_lvlind] = getLevelIndex(    comp_var,levels,    base_lvl);
        [comp_lat, comp_latind] = getLatitudeIndex( comp_var,YG,YC);
        [comp_lon, comp_lonind] = getLongitudeIndex(comp_var,XG,XC);
        
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
        
        base_var        = getVariableProperties(pinfo.vars);
        [level, lvlind] = getLevelIndex(    base_var,levels);
        [  lat, latind] = getLatitudeIndex( base_var,YG,YC);
        [  lon, lonind] = getLongitudeIndex(base_var,XG,XC);
        copyindices     = SetCopyID(fname);
        copy            = length(copyindices);
        
        pinfo.var_names      = base_var;
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
        var1                   = getVariableProperties(pinfo.vars);
        [var1_lvl, var1_lvlind] = getLevelIndex(    var1, levels);
        [var1_lat, var1_latind] = getLatitudeIndex( var1, YG, YC);
        [var1_lon, var1_lonind] = getLongitudeIndex(var1, XG, XC );
        
        disp('Getting information for the ''Y'' variable.')
        var2                   = getVariableProperties(pinfo.vars,    var1    );
        [var2_lvl, var2_lvlind] = getLevelIndex(    var2, levels, var1_lvl);
        [var2_lat, var2_latind] = getLatitudeIndex( var2, YG, YC, var1_lat);
        [var2_lon, var2_lonind] = getLongitudeIndex(var2, XG, XC, var1_lon);
        
        disp('Getting information for the ''Z'' variable.')
        var3                   = getVariableProperties(pinfo.vars,    var1    );
        [var3_lvl, var3_lvlind] = getLevelIndex(    var3, levels, var1_lvl);
        [var3_lat, var3_latind] = getLatitudeIndex( var3, YG, YC, var1_lat);
        [var3_lon, var3_lonind] = getLongitudeIndex(var3, XG, XC, var1_lon);
        
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

%% --------------------------------------------------------------------
%  All functions below this are intentionally scoped internal to main fcn.
%  --------------------------------------------------------------------

    function metadata = getVariableProperties(defvar)
        % With no arguments, getVariableProperties prompts the user to
        % specify a variable. The output is the variable and all the
        % coordinate arrays for that variable.
        % With one input argument, that argument is used as the default
        % when prompting. Otherwise, the first variable is usually the
        % default.
        
        if (nargin == 1)
            pgvar = defvar;
        else
            pgvar = pinfo.def_var;
        end
        
        % Build a string of potential variables to choose.
        str = sprintf(' %s ',pinfo.vars{1});
        for i = 2:pinfo.num_state_vars,
            str = sprintf('%s %s ',str,pinfo.vars{i});
        end
        
        fprintf('Default variable is ''%s'', if this is OK, <cr>;\n',pgvar)
        fprintf('If not, please enter one of: %s\n',str)
        varstring = input('(no syntax required)\n','s');

        if ~isempty(varstring), pgvar = varstring; end
        
        matcharray = strcmp(pgvar, pinfo.vars);
        if any(matcharray)
            pgvar = pinfo.vars{matcharray};
        else
            error('{%s} did not match anything in {%s}', pgvar, str)
        end
        
        % Must get all the metadata for this variable.
        % Read the 'original_coordinates' variable attribute, parse the coordinate
        % variable names, read the coordinate variables and stuff them into the
        % metadata structure. Easy.
        
        ncid = netcdf.open(pinfo.diagn_file,'NC_NOWRITE');

%        Variable Functions
%        ------------------
%        defVar           - Create NetCDF variable.
%        defVarChunking   - Set chunking layout.
%        defVarDeflate    - Set variable compression.
%        defVarFill       - Set fill parameters for variable.
%        defVarFletcher32 - Set checksum mode.
%        getVariableProperties           - Return data from NetCDF variable.
%        inqVar           - Return information about variable.
%        inqVarChunking   - Return chunking layout for variable.
%        inqVarDeflate    - Return variable compression information.
%        inqVarFill       - Return fill value setting for variable.
%        inqVarFletcher32 - Return checksum settings.
%        inqVarID         - Return ID associated with variable name.
%        putVar           - Write data to NetCDF variable.
%        renameVar        - Change name of NetCDF variable.
%        
%        Attribute Functions
%        -------------------
%        copyAtt          - Copy attribute to new location.
%        delAtt           - Delete NetCDF attribute.
%        getAtt           - Return NetCDF attribute.
%        inqAtt           - Return information about NetCDF attribute.
%        inqAttID         - Return ID of NetCDF attribute.
%        inqAttName       - Return name of NetCDF attribute.
%        putAtt           - Write NetCDF attribute.
%        renameAtt        - Change name of attribute.
 
       varid          = netcdf.inqVarID(ncid, pgvar);
       coordVarString = netcdf.getAtt(ncid, varid, 'original_coordinates');
       coordVars      = strsplit(coordVarString);
       metadata.varname = pgvar; % must initialize structure to something
 
       for i = 1:length(coordVars),
          varid    = netcdf.inqVarID(ncid, coordVars{i});
          
          % Could use original variable names as structure elements.
          % metadata.(coordVars{i}) = netcdf.getVariableProperties(ncid, varid);

          longName = netcdf.getAtt(ncid,varid,'long_name');
          switch lower(longName)
              case {'longitude'}
                  metadata.longitude = netcdf.getVar(ncid, varid);
              case {'latitude'}
                  metadata.latitude = netcdf.getVar(ncid, varid);
              case {'depth'}
                  metadata.level = netcdf.getVar(ncid, varid);
              case {'time'}
                  metadata.time = pinfo.time;
          end
       end
       
       netcdf.close(ncid);
       
       % TJH FIXME presumes nlon - nlat - nlev storage ...
       metadata.nlon = size(metadata.longitude,1);
       metadata.nlat = size(metadata.longitude,2);
       metadata.nlev = size(metadata.longitude,3);
       
    end % function getVariableProperties

%% --------------------------------------------------------------------

    function [time, timeind] = GetTime(pgvar, times, deftime)
        if (nargin == 3), time = deftime; else time = mean(times); end
        
        fprintf('Default time is %f, if this is OK, <cr>;\n',time)
        fprintf('If not, enter a time between %.4f and %.4f, we use the closest.\n', ...
            min(times),max(times))
        varstring = input('(no syntax required)\n','s');
        
        if ~isempty(varstring), time  = str2double(varstring); end
        
        d       = abs(time - times);    % crude distance
        ind     = find(min(d) == d);  % multiple minima possible
        timeind = ind(1);             % use the first one
        time    = times(timeind);
        
    end % GetTime

%% --------------------------------------------------------------------

    function lvlind = getLevelIndex(defaultLevel)
        if (nargin == 1)
            lvlind = defaultLevel;
        else
            lvlind = round(metadata.nlev/2);
        end

        fprintf('Default level (index) is  %d, if this is OK, <cr>;\n',lvlind)
        fprintf(' Enter an index between 1 and %d, inclusive ...\n', metadata.nlev)
        varstring = input('we''ll use the closest (no syntax required)\n','s');
        
        if ~isempty(varstring), lvlind = round(str2double(varstring)); end
        
    end % function getLevelIndex

%% --------------------------------------------------------------------

    function lonind = getLongitudeIndex(deflon)
        if (nargin == 1),
            lonind = deflon;
        else
            lonind = round(metadata.nlon/2);
        end
        
        fprintf('Default longitude (index) is %d, if this is OK, <cr>;\n',lonind)
        fprintf('If not, enter an index between 1 and %d, inclusive.\n', metadata.nlon)
        varstring = input('(no syntax required)\n','s');
        
        if ~isempty(varstring), lonind = round(str2double(varstring)); end
        
    end  % function getLongitudeIndex

%% --------------------------------------------------------------------

    function latind = getLatitudeIndex(deflat)
        if (nargin == 1),
            latind = deflat;
        else
            latind = round(metadata.nlat/2);
        end
        
        fprintf('Default latitude (index) is %d, if this is OK, <cr>;\n',latind)
        fprintf('If not, enter an index between 1 and %d, inclusive.\n',metadata.nlat)
        varstring = input('(no syntax required)\n','s');
        
        if ~isempty(varstring), latind = round(str2double(varstring)); end
        
    end  % function getLatitudeIndex

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

