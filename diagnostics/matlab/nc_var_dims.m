function vdims = nc_var_dims(fid,varname);
% Get the names of the coordinate variables for each 
% of the dimensions of the variable.
%
% vdims = nc_var_dims(fid,varname);
%
% fid        a netcdf file id (result of fid = netcdf(fname);
% varname    a variable names
% vdims      a cell array of the coordinate variables
%
% EXAMPLE:
% 
% fid   = netcdf('obs_seq.final.nc','nowrite');
% varname = 'RADIOSONDE_TEMPERATURE_guess';
% vdims = nc_var_dims(fid,varname);
% for i = 1:length(vdims)
%    disp(sprintf('variable %s dimension %d is %s',varname,i,vdims{i}))
% end

ncvobj = ncvar(varname,fid);
dims   = dim(ncvobj);
vdims  = [];

for j = 1:length(dims)
   vdims{j} = name(dims{j});
end
