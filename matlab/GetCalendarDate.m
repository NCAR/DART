function h = GetCalendarDate(seconds,days,calendartype)
%
% h = GetCalendarDate(seconds,days [,calendartype] )
%
% seconds, days are the DART times 
%
% EXAMPLE:
%
% mydate = GetCalendarDate(82761,148520);

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL:
% http://subversion.ucar.edu/DAReS/DART/trunk/diagnostics/matlab/fit_ens_mean_time.m
% $
% $Id$
% $Revision$
% $Date$

if (nargin < 2) 
   error('Must supply at least two arguments, seconds and days')
elseif (nargin ==2)
   calendartype = 'Gregorian';
end

switch lower(calendartype)
   case 'gregorian'
      mytime = datenum(1601,1,1) + days + seconds/86400;
      h = datestr(mytime);
      fprintf('DART time (%d s, %d d) is %s %s\n', ...
                     seconds,days,h,calendartype)
   case 'noleap'
      error('noleap not supported yet')
   case 'thirty_day_months'
      error('thirty_day_months not supported yet')
   case 'julian'
      error('julian not supported yet')
   case 'no_calendar'
      error('no_calendar not supported yet')
   case 'gregorian_mars'
      error('gregorian_mars not supported yet')
   otherwise
end

