function bob = get_varnames(fid)
% get_varnames returns JUST the 'atmospheric' variable names in the netCDF file - ie - not the coordinate variables.
%
% the result is a cell array of strings ... must use {} notation to address elements.
%
% EXAMPLE:
% fid = netcdf('obs_seq.final.nc','nowrite');
% varnames = get_varnames(fid);
% varnames{:}
% nvars = length(varnames);
% disp(sprintf('first atmospheric variable (of %d) is %s',nvars,varnames{1}))

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

% Use the netcdf primitives 'var','coord' to get complete list
% of (all) variables and the coordinate variables, respectively.

variables  =   var(fid);
coordvars  = coord(fid);

atmosvarinds = VariableIndices(coordvars,variables);

% coerce just the names into a cell array 

for i = 1:length(atmosvarinds)
   bob{i} = name(variables{atmosvarinds(i)});
end



function inds = VariableIndices(coordvars,variables)
% returns the indices of the atmospheric variables.

inds = [];

for i = 1:length(variables)

   varname    = name(variables{i});
   isatmosvar = 1;
   
   for j = 1:length(coordvars)
      if (strcmp( varname , name(coordvars{j}))), isatmosvar = 0; end 
      if (strcmp( varname ,      'time_bounds')), isatmosvar = 0; end 
      if (strcmp( varname ,     'region_names')), isatmosvar = 0; end 
      if (strcmp( varname ,     'CopyMetaData')), isatmosvar = 0; end 
   end

   if (isatmosvar > 0)
    inds = [inds i];
   end
end
