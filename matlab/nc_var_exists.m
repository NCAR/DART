function [variable_present, varid] = nc_var_exists(fname,varname)
%% NC_VAR_EXISTS Tests to see if the variable exists in a netCDF file - without throwing an error.
% If the variable exists, returns a 1, if it does not exist, it returns a 0
% Optionally returns the Matlab 'native' variable id ...
%
%
% EXAMPLE 1 (a failure):
% fname = 'obs_diag_output.nc';
% varname = 'variable_has_no_name';
% variable_present = nc_var_exists(fname,varname)
%
% EXAMPLE 2 (a success):
% fname = 'obs_diag_output.nc';
% varname = 'time';
% variable_present = nc_var_exists(fname,varname)
% if (variable_present), mytime = netcdf.getVar(ncid,varid); end
%
% EXAMPLE 3 (also a success):
% fname = 'obs_diag_output.nc';
% ncid  = netcdf.open(fname,'NOWRITE');
% varname = 'time';
% [variable_present, varid] = nc_var_exists(fname,varname);
% if (variable_present), mytime = netcdf.getVar(ncid,varid); end

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

finfo = ncinfo(fname);
variable_present = 0;
varid = [];

for ivar = 1:length(finfo.Variables)
   
   % TJH DEBUG  fprintf('Found Variable "%s" \n', vinfo.Variables(ivar).Name )
   if (strcmp(finfo.Variables(ivar).Name, deblank(varname)))
      % TJH DEBUG  fprintf('Matched Variable "%s" \n', vinfo.Variables(ivar).Name )
      variable_present = 1;
      ncid  = netcdf.open(fname,'NOWRITE');
      varid = netcdf.inqVarID(ncid,varname);
      netcdf.close(ncid);
      return
   end
end

% Purported Examples from Matlab
% vinfo       = ncinfo('example.nc','peaks');
% varSize     = vinfo.Size;
% disp(vinfo);
% hasUnLimDim = any([vinfo.Dimensions.Unlimited]);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
