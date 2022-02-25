function data = my_getVar(ncid,varid,FillValue)

% Usage data = my_getVar(ncid,Var_ind,FillValue);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

data = netcdf.getVar(ncid, varid);
if ( ~ isempty(FillValue) )
   data(data == FillValue) = NaN;
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
