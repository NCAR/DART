function dates = nc_read_time(fname,varid)
%% Attempts to interpret the time variable based on any existing 'units'
%  attribute. At present, the calendar is ignored.
%  Returns the time compatible with matlab's 'datenum' base.
%
% cdate = nc_read_time('example.nc','time');

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

dates = [];
timeorigin = 0;

times     = ncread(fname,varid);
timeunits = nc_read_att(fname,varid,'units');
if (~ isempty(timeunits))
    timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
    if (timebase(1) > 0000)
        timeorigin = datenum(timebase(1),timebase(2),timebase(3));
    end
end
dates      = times + timeorigin;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
