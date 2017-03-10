function pinfo = GetMPAS_ATMInfo(pinfo_in,fname,routine)
%% GetMPAS_ATMInfo   prepares a structure of information needed by the subsequent "routine"
%               The information is gathered via rudimentary "input" routines.
%               This function is NOT meant to be called directly by the user.
%
% pinfo = GetMPAS_ATMInfo(pinfo_in, fname, routine);
%
% pinfo_in  Name of existing pinfo struct, e.g. output from CheckModelCompatibility
% fname     Name of the DART netcdf file - usually [Prior,Posterior]_Diag.nc
% routine   name of subsequent plot routine.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

pinfo = pinfo_in;
model = ncreadatt(fname, '/', 'model');

if strcmpi(model,'mpas_atm') ~= 1
   error('Not so fast, this is not a MPAS_ATM model.')
end
if strcmpi(pinfo_in.model,model) ~= 1
   error('Not so fast, this is not a MPAS_ATM plotting structure.')
end

%% Get the domain information.

varexist(fname, {'lonCell','latCell','areaCell'})

pinfo.area     = ncread(fname,'areaCell');
pinfo.lonCell  = ncread(fname,'lonCell');
pinfo.latCell  = ncread(fname,'latCell');
pinfo.lonunits = ncreadatt(fname,'lonCell','units');
pinfo.latunits = ncreadatt(fname,'latCell','units');

%% code for each plot type

switch lower(deblank(routine))

    case {'plottotalerr'}

      pinfo.fname      = fname;

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar            = GetVarString(pinfo_in.vars);
      cellind          = GetClosestCell(pgvar, pinfo.latCell, pinfo.lonCell);
      [level, lvlind]  = GetLevel(fname, pgvar);

      pinfo.fname      = fname;
      pinfo.var        = pgvar;
      pinfo.level      = level;
      pinfo.levelindex = lvlind;
      pinfo.cellindex  = cellind;
      pinfo.longitude  = pinfo.lonCell(cellind);
      pinfo.latitude   = pinfo.latCell(cellind);

   case 'plotcorrel'

      disp('Getting information for the ''base'' variable.')
      base_var                 = GetVarString(pinfo_in.vars);
      base_cellind             = GetClosestCell(base_var, pinfo.latCell, pinfo.lonCell);
      [base_time, base_tmeind] = GetTime(pinfo.time);
      [base_lvl,  base_lvlind] = GetLevel(fname, base_var);

      disp('Getting information for the ''comparison'' variable.')
      comp_var                 = GetVarString(pinfo_in.vars, base_var);
      [comp_lvl, comp_lvlind]  = GetLevel(fname, comp_var, base_lvl);

      pinfo.fname          = fname;
      pinfo.base_var       = base_var;
      pinfo.comp_var       = comp_var;
      pinfo.base_time      = base_time;
      pinfo.base_tmeind    = base_tmeind;
      pinfo.base_lvl       = base_lvl;
      pinfo.base_lvlind    = base_lvlind;
      pinfo.base_cellindex = base_cellind;
      pinfo.base_lat       = pinfo.latCell(base_cellind);
      pinfo.base_lon       = pinfo.lonCell(base_cellind);
      pinfo.comp_lvl       = comp_lvl;
      pinfo.comp_lvlind    = comp_lvlind;

   case 'plotvarvarcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVarString(pinfo_in.vars);
       base_cellind            = GetClosestCell(base_var, pinfo.latCell, pinfo.lonCell);
      [base_time, base_tmeind] = GetTime(pinfo.time);
      [base_lvl , base_lvlind] = GetLevel(fname, base_var);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVarString(pinfo_in.vars,      base_var);
       comp_cellind           = GetClosestCell(comp_var, pinfo.latCell, pinfo.lonCell, base_cellind);
      [comp_lvl, comp_lvlind] = GetLevel(fname, comp_var, base_lvl);

      pinfo.fname          = fname;
      pinfo.base_var       = base_var;
      pinfo.comp_var       = comp_var;
      pinfo.base_time      = base_time;
      pinfo.base_tmeind    = base_tmeind;
      pinfo.base_lvl       = base_lvl;
      pinfo.base_lvlind    = base_lvlind;
      pinfo.base_cellindex = base_cellind;
      pinfo.base_lat       = pinfo.latCell(base_cellind);
      pinfo.base_lon       = pinfo.lonCell(base_cellind);
      pinfo.comp_lvl       = comp_lvl;
      pinfo.comp_lvlind    = comp_lvlind;
      pinfo.comp_cellindex = comp_cellind;
      pinfo.comp_lat       = pinfo.latCell(comp_cellind);
      pinfo.comp_lon       = pinfo.lonCell(comp_cellind);

   case 'plotsawtooth'

       pgvar                      = GetVarString(pinfo_in.vars);
      [level, lvlind]             = GetLevel(pinfo_in.prior_file, pgvar);
       cellind                    = GetClosestCell(pgvar, pinfo.latCell, pinfo.lonCell);
      [copyindices, copymetadata] = SetCopyID2(pinfo_in.prior_file);
      copy                        = length(copyindices);

      pinfo.truth_file     = pinfo_in.truth_file;
      pinfo.prior_file     = pinfo_in.prior_file;
      pinfo.posterior_file = pinfo_in.posterior_file;
      pinfo.var_names      = pgvar;
      pinfo.level          = level;
      pinfo.levelindex     = lvlind;
      pinfo.cellindex      = cellind;
      pinfo.latitude       = pinfo.latCell(cellind);
      pinfo.longitude      = pinfo.lonCell(cellind);
      pinfo.copies         = copy;
      pinfo.copyindices    = copyindices;

      % Hui has a habit of cutting out the rest of the ensemble members
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
      [var1_lvl, var1_lvlind] = GetLevel(fname,var1);
       var1_cellind           = GetClosestCell(var1, pinfo.latCell, pinfo.lonCell);

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVarString(pinfo_in.vars,        var1    );
      [var2_lvl, var2_lvlind] = GetLevel(fname,var2,var1_lvl);
       var2_cellind           = GetClosestCell(var2, pinfo.latCell, pinfo.lonCell, var1_cellind);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVarString(pinfo_in.vars,        var1    );
      [var3_lvl, var3_lvlind] = GetLevel(fname,var3,var1_lvl);
       var3_cellind           = GetClosestCell(var3, pinfo.latCell, pinfo.lonCell, var1_cellind);

      % query for ensemble member string
      metastrings = ncread(fname,'CopyMetaData');
      if(size(metastrings,2) == 1), metastrings = metastrings'; end
      metadata    = cellstr(metastrings);
      ens_mem     = strtrim(metadata{1});
      str1 = sprintf('Input ensemble member metadata STRING. <cr> for ''%s''\n',ens_mem);
      s1   = input(str1,'s');
      if ~ isempty(s1), ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''\n','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

      pinfo.fname          = fname;
      pinfo.var1name       = var1;
      pinfo.var2name       = var2;
      pinfo.var3name       = var3;
      pinfo.var1_lvl       = var1_lvl;
      pinfo.var1_lvlind    = var1_lvlind;
      pinfo.var1_cellindex = var1_cellind;
      pinfo.var1_lat       = pinfo.latCell(var1_cellind);
      pinfo.var1_lon       = pinfo.lonCell(var1_cellind);
      pinfo.var2_lvl       = var2_lvl;
      pinfo.var2_lvlind    = var2_lvlind;
      pinfo.var2_cellindex = var2_cellind;
      pinfo.var2_lat       = pinfo.latCell(var2_cellind);
      pinfo.var2_lon       = pinfo.lonCell(var2_cellind);
      pinfo.var3_lvl       = var3_lvl;
      pinfo.var3_lvlind    = var3_lvlind;
      pinfo.var3_cellindex = var3_cellind;
      pinfo.var3_lat       = pinfo.latCell(var3_cellind);
      pinfo.var3_lon       = pinfo.lonCell(var3_cellind);
      pinfo.ens_mem        = ens_mem;
      pinfo.ltype          = ltype;

    otherwise

        error('%s unsupported for %s',routine, pinfo.model)

end



function pgvar = GetVarString(prognostic_vars, defvar)
%% --------------------------------------------------------------------
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
%% --------------------------------------------------------------------
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



function [level, lvlind] = GetLevel(fname, pgvar, deflevel)
%% --------------------------------------------------------------------
%

varinfo  = nc_getvarinfo(fname, pgvar);
vertdim  = [];

for idim = 1:length(varinfo.Dimension),
   if (strncmpi(varinfo.Dimension{idim},'nvert',5) )
      vertdim = idim;
      break
   end
end

if (isempty(vertdim))

   fprintf('''%s'' has no vertical level, imply using 1.\n',pgvar)
   level  = 1;
   lvlind = 1;

else

   nlevels = varinfo.Size(vertdim);
   levels  = 1:nlevels;
   if (nargin == 3)
      level = deflevel;
   else
      level = round(nlevels/2);
   end

   fprintf('Default level index is # %d, if this is OK, <cr>;\n',level)
   fprintf('If not, enter a level index between 1 and %d, inclusive ...\n',nlevels)
   varstring = input('we''ll use the closest (no syntax required)\n','s');

   if ~isempty(varstring), level = str2num(varstring); end

   d      = abs(level - levels);  % crude distance
   ind    = find(min(d) == d);    % multiple minima possible
   lvlind = ind(1);               % use the first one
   level  = levels(lvlind);

end



function cellindex = GetClosestCell(pgvar, lats, lons, defcellind)
%% --------------------------------------------------------------------
%  The MPAS_ATM domain is an unstructured grid.
%  Could/Should see if the variable is a Cell-centric quantity.

ncells = length(lats);

% Determine a sensible default.

if (nargin > 3)
   cellind = defcellind;
else
   cellind = round(ncells/2);
end
lat = lats(cellind);
lon = lons(cellind);

% Ask if they like the default.

fprintf('There are %d cells for %s\n',ncells, pgvar)
fprintf('Default cell index is %d ... lat/lon(%f,%f)\n',cellind,lat,lon)
fprintf('If this is OK, <cr>; if not, enter desired lat/lon pair ...\n')
varstring = input('we''ll use the closest (no syntax required)\n','s');

if ~isempty(varstring)
   nums = str2num(varstring);
   if (length(nums) ~= 2)
      error('Did not get two indices for the lat lon pair.')
   end
   lat = nums(1);
   lon = nums(2);
end

% For simplicity I will pick the grid cell whose center is closest
% to the point of interest. This is not a perfect solution, should really
% check the neighboring cells to see if it is INSIDE the cell,
% but if you are interested in that level of accuracy, you're on your own.

londiff   = lons - lon;
latdiff   = lats - lat;
dist      = sqrt(londiff.^2 + latdiff.^2);
cellindex = find(min(dist(:)) == dist);



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
