function slab = get_hyperslab(varargin)
%% DART:get_hyperslab gets an arbitrarily-shaped variable from a netCDF diagnostic file.
% Retrieves a particular copy of a state vector from a file whose
% full or relative path is specified in the file argument.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

for i = 1:2:nargin,
   eval(sprintf('pinfo.%s = varargin{i+1};',varargin{i}))
end

if ( exist(pinfo.fname,'file') ~= 2 ), error('%s does not exist.',pinfo.fname); end

[start, count] = GetNCindices(pinfo,'fname',pinfo.varname);
slab           = ncread(pinfo.fname, pinfo.varname, start, count);

if (sum(isfinite(slab(:))) == 0)
   pinfo
   error('%s %s has all missing values ... exiting.', pinfo.fname, pinfo.varname)
end

% The mere presence of the field is enough
if (isfield(pinfo,'squeeze'))
    slab = squeeze(slab);
end

% The mere presence of the field is enough
% permute is the equivalent of an N-D transpose
if (isfield(pinfo,'permute'))
    rank = length(size(slab));
    slab = permute(slab,rank:-1:1);
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
