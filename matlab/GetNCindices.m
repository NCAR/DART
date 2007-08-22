function [corner, endpnt] = GetNCindices(pinfo, whichfile, varname);
% GETNCindices returns a corner,endpnt array for use with getnc.
% At present, all times, all copies for a specific level,lat,lon.
% Does not assume anything about the dimension of the variable.
%
% USAGE:
% [corner, endpnt] = GetNCindices(pinfo, whichfile, varname);
% whichfile          is a character string specifying which 
%                    component of 'pinfo' will be used.
%                    ['prior','posterior','truth']
% varname            is the netcdf variable being extracted.

% The structure 'pinfo' must have the following components:
%pinfo.latindex
%pinfo.lonindex
%pinfo.levelindex
% and the value  

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

% GetNCindices replaces the following hardwired piece of code.
%
%      if ( strcmp(lower(vname),'ps') ==1 ) %  PS(time, copy, lat, lon)
%         corner = [  1 NaN pinfo.latindex pinfo.lonindex];
%         endpnt = [ -1 NaN pinfo.latindex pinfo.lonindex];
%      else % U(time, copy, lev, lat, lon)
%         corner = [  1 NaN pinfo.levelindex pinfo.latindex pinfo.lonindex ];
%         endpnt = [ -1 NaN pinfo.levelindex pinfo.latindex pinfo.lonindex ];
%      end

switch lower(whichfile)
   case 'prior'
      fname = pinfo.prior_file;
   case 'posterior'
      fname = pinfo.posterior_file;
   case 'truth'
      fname = pinfo.truth_file;
   case 'diagn'
      fname = pinfo.diagn_file;
   otherwise
      fname = pinfo.fname;
end

corner = [];
endpnt = [];

% This block uses the lowest level netcdf operator I know of.

[cdfid, rcode] = ncmex('OPEN',fname,'NC_NOWRITE');
[varid, rcode] = ncmex('VARID', cdfid, varname);
[var, datatype, ndims, dim, natts, status] = ncmex('VARINQ', cdfid, varid);

% var      is the name of the variable as it is in the netCDF file
% datatype is 
% ndims    is the 'rank' of the variable 2D, 3D, 4D, etc.
% dims     is an array of length 'ndims' specifying the 
%          dimension ID's used - IN ORDER.
% natts    is the number of attributes for the variable
% status   is obvious

% make sure it has a time and copy index

% loop over all of the variables dimensions and 
% build up the corner/endpoint arrays

for i = 1:ndims 

   [dimname, length, status] = ncmex('DIMINQ', cdfid, dim(i));

   % the coordinate variable is the same name as the dimension
   % some grids have multiple vertical levels so there is no one
   % predictable coordinate variable name. By convention, these
   % variables should be tagged with a 'cartesian_axis' attribute.
   % ditto for lat, lon ... (on staggered grids, for example)

%   [var, datatype, nddims, ddim, ndatts, status] = ncmex('VARINQ', cdfid, dimname);
%   [datatype, len, status] = ncmex('ATTINQ', cdfid, var, 'cartesian_axis');
   
   [datatype, len, status] = ncmex('ATTINQ', cdfid, dimname, 'cartesian_axis');
   
   if (status > 0) 
      % There is a 'cartesian_axis' attribute, and it should be one of:
      % [X,Y,Z,T]
      [value, stat2] = ncmex('ATTGET', cdfid, dimname, 'cartesian_axis');
      dimname = value;
   else
      % Then there is no 'cartesian_axis' attribute and the best we can
      % hope for is a standard name [time,copy,lat,lon,lev]
   end

   switch lower(dimname)
      case {'time','t'}
         corner = [corner  1];
         endpnt = [endpnt length];
      case 'copy'
         corner = [corner  1];
         endpnt = [endpnt length];
      case {'lat','y','tmpj'}
         corner = [corner pinfo.latindex ];
         endpnt = [endpnt pinfo.latindex ];
      case {'lon','x','tmpi'}
         corner = [corner pinfo.lonindex ];
         endpnt = [endpnt pinfo.lonindex ];
      case {'lev','z'}
         corner = [corner pinfo.levelindex ];
         endpnt = [endpnt pinfo.levelindex ];
   end

end

ncmex('CLOSE',fname);

