function timestring = set_time_axis(whichone,bincenters,DateForm)
%%
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if strncmpi(whichone,'x',1)
    ax = 'x';
else
    ax = 'y';
end

ndays = max(bincenters(:)) - min(bincenters(:)) + 1;

if (bincenters(1) > 1000) && (ndays > 5)
    if  strncmpi(DateForm,'default',7)
        dateform = 'mm/dd';
    else
        dateform = DateForm;
    end
    datetick('x',dateform,'keeplimits','keepticks');
    monstr = datestr(bincenters(1),21);
    timestring = sprintf('%s start',monstr);
elseif (bincenters(1) > 1000)
    if  strncmpi(DateForm,'default',7)
        dateform = 'dd HH:MM';
    else
        dateform = DateForm;
    end
    datetick(ax,dateform,'keeplimits')
    monstr = datestr(bincenters(1),21);
    timestring = sprintf('%s start',monstr);
else
    timestring = 'days';
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

