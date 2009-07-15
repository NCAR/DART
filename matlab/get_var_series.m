function var_vec = get_var_series(fname, varname, copyindex, state_var, tstart, tcount)
%GET_VAR_SERIES Gets a particular copy of a state variable from netcdf file
%
% Retrieves a particular copy of a state variable from a file whose
% full or relative path is specified in the file argument.
%
% Example 1:
% fname     = '../work/Prior_Diag.nc';
% varname   = 'state';      % State Variable
% copyindex = 8;            % copy index (Ensemble Member)
% state_var = 3;            % state variable index
% var_vec   = get_var_series(fname, varname, copyindex, state_var);

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

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

if (nargin == 4)
  tstart =  1;
  tcount = -1;
end

num_times   = dim_length(fname,'time');
num_copies  = dim_length(fname,'copy');
num_vars    = dim_length(fname,'StateVariable');

if (copyindex > num_copies ) 
    fprintf('%s only has %d ''copies/Ensemble members''\n',fname,num_copies)
    error('you wanted copy %d ', copyindex)
end
if (state_var > num_vars) 
   fprintf('%s only has %d %s variables\n',fname,num_vars,varname)
   error('you wanted variable %d ', state_var)
end

% Get only the appropriate copy of the state and return

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.stateindex = state_var;
[start, count]    = GetNCindices(myinfo, 'diagn', varname);

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstart-1;
         count(i) = tcount;
         break
      otherwise
   end
end

var_vec = nc_varget(fname, varname, start, count);

if (sum(isfinite(var_vec(:))) == 0) 
   error('%s %s copy %d index %d has all missing values ... exiting.', ...
        fname,varname,copyindex,state_var)
end


function x = dim_length(fname,dimname)

bob = nc_getdiminfo(fname,dimname);
x   = bob.Length;

