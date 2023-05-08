function [] = write_netcdf_att(ncname,varname,fillval,units, ...
    long,standard,internal,cell)
% function write_netcdf_att: loads all years of tower met forcing (PLUMBER2 format)
% Input:  'ncname' netcdf file name and path
%         'varname' is the CAM met name
%         'fillval' FillValue
%         'units'
%         'long,standard,internal' all attribute naming conventions
%         'cell': optional argument
% Output: N/A


ncwriteatt(ncname,varname,'_FillValue',fillval)
ncwriteatt(ncname,varname,'units',units)
ncwriteatt(ncname,varname,'long_name',long)
ncwriteatt(ncname,varname,'standard_name',standard)
ncwriteatt(ncname,varname,'internal_dname',internal)
if exist('cell','var')
    ncwriteatt(ncname,varname,'cell_methods',cell)
end

end
