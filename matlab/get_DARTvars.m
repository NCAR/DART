function bob = get_DARTvars(fname)
%% get_DARTvars returns just the variable names in the netCDF file that are likely to
%be DART variables.
%
% the result is a cell array of strings ... must use {} notation to address elements.
%
% EXAMPLE:
% fname    = 'obs_seq.final.nc';
% DARTvars = get_DARTvars(fname);
% DARTvars{:}
% nvars = length(DARTvars);
% disp(sprintf('first atmospheric variable (of %d) is %s',nvars,DARTvars{1}))

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fileinfo  = nc_info(fname);
nvars     = length(fileinfo.Dataset);
isDARTvar = zeros(nvars,1);

for i = 1:nvars

   % Anything with a 'copy' dimension is probably a DART state vector variable.

   dimnames = fileinfo.Dataset(i).Dimension;

   if (any(strcmp(dimnames,'copy'))), isDARTvar(i) = 1; end

   % Reject the obvious coordinate variables and some metadata ones

   varname    = fileinfo.Dataset(i).Name;
   if (         nc_iscoordvar(fname,varname)), isDARTvar(i) = 0; end
   if (strcmp( varname ,     'CopyMetaData')), isDARTvar(i) = 0; end

end

if (sum(isDARTvar) == 0)
   error('No DART state variables in %s',fname)
end

% coerce just the names into a cell array

varind = 0;
for i = 1:nvars
   if (isDARTvar(i) > 0)
      varind = varind + 1;
      bob{varind} = fileinfo.Dataset(i).Name;
   end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
