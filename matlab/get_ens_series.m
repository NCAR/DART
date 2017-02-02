function ens = get_ens_series(fname, varname, state_var_index, tstartind, tendind)
%% GET_ENS_SERIES: Returns matrix of time series for all members of ensemble for a variable
%
% the rows of the matrix correspond to time,
% the columns of the matrix correspond to ensemble members
%
% fname = 'Prior_Diag.nc';
% varname = 'state';
% state_var_index = 3;
% ens = get_ens_series(fname,varname,state_var_index);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

disp('get_ens_series() is deprecated, use get_hyperslab() instead.')

if (nargin == 3) 
  tstartind =  1;
  diminfo   = nc_getdiminfo(fname,'time');
  tendind   = diminfo.Length;
end

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

model      = nc_attget(fname,nc_global,'model');
varinfo    = nc_getvarinfo(fname,varname);
num_copies = varinfo.Size(2);
num_vars   = varinfo.Size(3);

if ( ~ strcmp( varinfo.Dimension{1}, 'time') )
   fprintf('%s first dimension ( %s ) is not ''time''\n',fname,varinfo.Dimension{1})
end
if ( ~ strcmp( varinfo.Dimension{2}, 'copy') )
   fprintf('%s second dimension ( %s ) is not ''copy''\n',fname,varinfo.Dimension{2})
end
if (state_var_index > num_vars)
   fprintf('%s only has %d %s variables\n',fname,num_vars,varname)
   fprintf('you wanted variable %d\n', state_var_index)
end

[ens_num, copyindices] = get_ensemble_indices(fname);

% Get the whole thing and then return the ones we want.
% This is usually not too bad, as there are usually many more
% ensemble members than "mean" and "spread" (the two members
% we are NOT interested in for this function).

myinfo.diagn_file = fname;
myinfo.stateindex = state_var_index;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

state_vec = nc_varget(fname,varname,start,count);
ens       = state_vec(:,copyindices);

fprintf('Read %d ensemble members for variable %d in %s\n', ...
             ens_num, state_var_index,fname);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
