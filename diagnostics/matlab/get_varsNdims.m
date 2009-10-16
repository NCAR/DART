function [y, ydims] = get_varsNdims(fname);
% Get the dimension (strings) for each atmospheric variable.
% [y, ydims] = get_vars_dims(fname);
%
% fname     a netcdf file name
%
% y       a cell array of variable names
% ydims   a cell array of the concatenated dimension names 
%
% EXAMPLE:
% 
% fname      = 'obs_seq.final.nc';
% [y, ydims] = get_varsNdims(fname);
%
% >> plotdat.allvarnames{20}  
%
%    AIRCRAFT_U_WIND_COMPONENT_guess
%
% >> plotdat.allvardims{20}
%    region plevel copy time

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

ALLvarnames = get_varnames(fname);
Nvarnames   = length(ALLvarnames);

for i = 1:Nvarnames

   varname = ALLvarnames{i};
   varinfo = nc_getvarinfo(fname,varname);

   y{i}     = varname;
   ydims{i} = sprintf('%s ',varinfo.Dimension{:});

end
