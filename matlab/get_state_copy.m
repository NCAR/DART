function state_vec = get_state_copy(fname, varname, copyindex, tstartind, tcount)
%% GET_STATE_COPY  Gets a particular copy (one ensemble member) of state from netcdf file
% Retrieves a particular copy of a state vector from a file whose
% full or relative path is specified in the file argument.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

disp('get_state_copy() is deprecated, use get_hyperslab() instead.')


if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

if (nargin == 3)
  tstartind =  1;
  diminfo   = nc_getdiminfo(fname,'time');
  tcount    = diminfo.Length;
end

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
[start, count]    = GetNCindices(myinfo, 'diagn', varname);

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tcount;
         break
      otherwise
   end
end

state_vec = nc_varget(fname, varname, start, count);

if (sum(isfinite(state_vec(:))) == 0)
   error('%s %s copy %d has all missing values ... exiting.', ...
        fname, varname, copyindex )
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
