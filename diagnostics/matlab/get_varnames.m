function bob = get_varnames(fname)
%% get_varnames returns JUST the 'atmospheric' variable names in the netCDF file - ie - not the coordinate variables.
%
% the result is a cell array of strings ... must use {} notation to address elements.
%
% EXAMPLE:
% fname    = 'obs_seq.final.nc';
% varnames = get_varnames(fname);
% varnames{:}
% nvars = length(varnames);
% disp(sprintf('first atmospheric variable (of %d) is %s',nvars,varnames{1}))

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fileinfo    = nc_info(fname);
Nvarnames   = length(fileinfo.Dataset);

inds = [];

for i = 1:Nvarnames

   varname    = fileinfo.Dataset(i).Name;
   isatmosvar = 1;

   % Reject the obvious coordinate variables and some that are
   % specific to DART

   if (         nc_iscoordvar(fname,varname)), isatmosvar = 0; end 
   if (strcmp( varname ,      'time_bounds')), isatmosvar = 0; end 
   if (strcmp( varname ,     'region_names')), isatmosvar = 0; end 
   if (strcmp( varname ,     'CopyMetaData')), isatmosvar = 0; end 
   if (strcmp( varname , 'ObservationTypes')), isatmosvar = 0; end

   % keep track of the 'good' variables
   if (isatmosvar > 0)
      inds = [inds i];
   end
end

if (isempty(inds))
   error('No atmospheric variables in %s',fname)
end

% coerce just the names into a cell array 

for i = 1:length(inds)
   bob{i} = fileinfo.Dataset(inds(i)).Name;
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

