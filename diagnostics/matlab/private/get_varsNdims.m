function [y, ydims] = get_varsNdims(fname)
%% Get the dimension (strings) for each atmospheric variable.
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
% >> y{20}  
%
%    AIRCRAFT_U_WIND_COMPONENT_guess
%
% >> ydims{20}
%    region plevel copy time

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

ALLvarnames = get_varnames(fname);
Nvarnames   = length(ALLvarnames);

y     = cell(Nvarnames,1);
ydims = cell(Nvarnames,1);

for i = 1:Nvarnames

   varname = ALLvarnames{i};
   y{i}    = varname;
   
   varinfo = ncinfo(fname,varname);
   ydims{i} = varinfo.Dimensions(1).Name;
   
   rank = length(varinfo.Size);
    
   for idim = 2:rank
      ydims{i} = sprintf('%s %s',ydims{i}, varinfo.Dimensions(idim).Name);
   end

end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
