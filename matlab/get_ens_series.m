function ens = get_ens_series(fname, varname, state_var_index, tstart, tend)
%GET_ENS_SERIES: Returns matrix of time series for all members of ensemble for a variable
%
% the rows of the matrix correspond to time,
% the columns of the matrix correspond to ensemble members
%
% fname = 'Prior_Diag.nc';
% varname = 'state';
% state_var_index = 3;
% ens = get_ens_series(fname,varname,state_var_index);

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

if (nargin == 3) 
  tstart = -1;
  tend = -1;
end

f = netcdf(fname);
model      = f.model(:);
var_atts   = dim(f{varname});       % cell array of dimensions for the var
num_copies = length(var_atts{2});
num_vars   = length(var_atts{3});

if ( ~ strcmp( name(var_atts{1}), 'time') )
    disp( sprintf('%s first dimension ( %s ) is not ''time''',fname,name(var_atts{1})))
end
if ( ~ strcmp( name(var_atts{2}), 'copy') )
    disp( sprintf('%s second dimension ( %s ) is not ''copy''',fname,name(var_atts{2})))
end
if (state_var_index > num_vars)
   disp( sprintf('%s only has %d %s variables',fname,num_vars,varname))
   error(sprintf('you wanted variable %d ', state_var_index))
end
close(f);

metadata    = getnc(fname,'CopyMetaData');           % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's
if ( isempty(copyindices) )
   disp(sprintf('%s has no valid ensemble members',fname))
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   disp(sprintf('%s claims to have %d copies',fname, num_copies))
   error('netcdf file has no ensemble members.')
end
ens_num     = length(copyindices);

% Get the whole thing and then return the ones we want.
% This is usually not too bad, as there are usually many more
% ensemble members than "mean" and "spread" (the two members
% we are NOT interested in for this function).

state_vec = getnc(fname,varname, [tstart, -1, state_var_index], ...
                                 [tend,   -1, state_var_index]);
% getnc always squeezes out the singleton last dimension.
ens       = state_vec(:,copyindices);

disp(sprintf('Read %d ensemble members for variable %d in %s', ...
             ens_num, state_var_index,fname));

