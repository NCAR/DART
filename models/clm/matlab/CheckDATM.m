%%
%
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

dirname = '/glade/proj3/DART/CAM_DATM';

% for iyear = [1999 2008 2009 2010]
for iyear = 2009:2009
for imem = 1:80

   filename = sprintf('%s/%d/CAM_DATM-%02d.cpl.ha2x1dx6h.%d.nc', ...
                      dirname,iyear,imem,iyear);

   if ( exist(filename,'file') ~= 2)
      error('%s does not exist',filename)
   end

   times       = nc_varget(filename,'time');
   timeunits   = nc_attget(filename,'time','units');
   timebase    = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
   timeorigin  = datenum(timebase(1),timebase(2),timebase(3));
   ttimes      = times + timeorigin;

   ntimes = length(ttimes);
   deltat = diff(ttimes);
   delta  = median(deltat);
   inds   = find(deltat > delta);

   if ((length(ttimes) == 365*4) && isempty(inds)) 
      fprintf('%s has all the timesteps.\n',filename)
   else
      fprintf('WARNING: %s has %d timesteps, not %d.\n',filename,ntimes,365*4)
      disp('WARNING: problem occurrs around')
      datestr(ttimes(inds))
   end

end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
