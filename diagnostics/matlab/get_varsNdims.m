function [y, ydims] = get_varsNdims(fid);
% Get the dimension (strings) for each atmospheric variable.
% [y, ydims] = get_vars_dims(fid);
%
% fid     a netcdf file id (result of fid = netcdf(fname);
%
% y       a cell array of variable names
% ydims   a cell array of the (appended) dimension names 
%
% EXAMPLE:
% 
% fid = netcdf('obs_seq.final.nc','nowrite');
% [y, ydims] = get_varsNdims(fid);
% nvars = length(y);
% disp(sprintf('variable %s has named dimensions %s',y{1},ydims{1}))

ALLvarnames = get_varnames(fid);
Nvarnames   = length(ALLvarnames);

ix = 0;

for i = 1:Nvarnames

   ncvobj = ncvar(ALLvarnames{i},fid);
   dims   = dim(ncvobj);              % cell array ... 
   dimnames = [];
   for j = 1:length(dims)
      dimnames = [dimnames ' ' name(dims{j})];
   end

   y{i}     = ALLvarnames{i};
   ydims{i} = dimnames;

end
