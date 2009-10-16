function vdims = nc_var_dims(ncfname,varname);
% Get the names of the coordinate variables for each 
% of the dimensions of the variable.
%
% vdims = nc_var_dims(ncfname,varname);
%
% ncfname    file name of a netcdf file
% varname    a variable names
% vdims      a cell array of the coordinate variables
%
% EXAMPLE:
% 
% ncfname = 'obs_seq.final.nc';
% varname = 'RADIOSONDE_TEMPERATURE_guess';
% vdims   = nc_var_dims(ncfname,varname);
% for i = 1:length(vdims)
%    disp(sprintf('variable %s dimension %d is %s',varname,i,vdims{i}))
% end

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

varinfo = nc_getvarinfo(ncfname,varname);
vdims   = varinfo.Dimension;
