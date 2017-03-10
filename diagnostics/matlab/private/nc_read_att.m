function value = nc_read_att(fname,varid,attname)
%% If the attribute exists, return the value; if not return an empty object.
%
% This differs from Matlab's intrinsic ncreadatt() in that it does not error
% out if the attribute does not exist.
%
% Some examples:
%
% cdate = nc_read_att('example.nc','/','creation_date');
% units = nc_read_att('example.nc','temperature','units');
% uhoh  = nc_read_att('example.nc','temperature','something_that_does_not_exit');
% if isempty(uhoh), fprintf('no such attribute but life goes on.\n'); end

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

value = [];

if (varid == '/')
    vinfo = ncinfo(fname);
else
    vinfo = ncinfo(fname,varid);
end

for iatt = 1:length(vinfo.Attributes)
    if (strcmp(vinfo.Attributes(iatt).Name, deblank(attname)))
        value = vinfo.Attributes(iatt).Value;
        return
    end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
