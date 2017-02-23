function [dimnames, dimsizes] = nc_var_dims(ncfname,varname)
%% Get the names of the coordinate variables for each 
% of the dimensions of the variable.
%
% dimnames = nc_var_dims(ncfname,varname);
%
% ncfname    file name of a netcdf file
% varname    a variable names
% dimnames      a cell array of the coordinate variables
%
% EXAMPLE:
% 
% ncfname = 'obs_seq.final.nc';
% varname = 'RADIOSONDE_TEMPERATURE_guess';
% dimnames   = nc_var_dims(ncfname,varname);
% for i = 1:length(dimnames)
%    disp(sprintf('variable %s dimension %d is %s',varname,i,dimnames{i}))
% end

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

dimnames = [];
dimsizes = [];

[present, ~] = nc_var_exists(ncfname,varname);

if ( present )
    varinfo  = ncinfo(ncfname,varname);
    dimsizes = varinfo.Size;
    rank     = length(dimsizes);
    dimnames = cell(1,rank);
    for idim = 1:rank
        dimnames{idim} = varinfo.Dimensions(idim).Name;
    end
else
    error('%s does not have a %s variable.',ncfname, varname)
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
