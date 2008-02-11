function h = SimpleMap(fname,var,time,copystring,level)
% SimpleMap plots a horizontal slice from a DART netcdf file.
% Right now, the hard assumption is the netcdf variable has
% order [time, copy, lat, lon (,lev)].
%
% Since WRF does not have tidy coordinate variables for lat/lon,
% this script will not directly support WRF/DART netCDF files.
%
% EXAMPLE:
%
% fname  = 'PriorDiag.nc';
% var    = 'T';      % 'PS','T','US','VS','Q','CLDLIQ','CLDICE' ...
% time   = 1;        % index into time array
% copy   = 'ensemble mean';  % anything from ncdump -v CopyMetaData
% level  = 1;        % index into level array
% h = SimpleMap(fname,var,time,copy,level);

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
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

% Get some information from the netCDF file.

plotdat = netCDFinfo(fname,var);

% Figure out the starting date

timebase   = sscanf(plotdat.timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
timearr    = getnc(fname,'time');
timestring = datestr(timearr(1) + timeorigin);

% create the hyperslab indices ... 
% most are [time,copy,lat,lon,lev] but some (i.e. PS)
% have no 'lev' dimension.
% There is a HARD assumption about variable shape here ...

copyindex = get_copy_index(fname,copystring);

if ( length(plotdat.dimcells) < 5 )
   corner    = [time,copyindex,-1,-1];
   end_point = [time,copyindex,-1,-1];
else 
   corner    = [time,copyindex,-1,-1,level];
   end_point = [time,copyindex,-1,-1,level];
end
   
datmat = getnc(fname,var,corner,end_point);

% Create the plot.
% The figure window is partitioned into two parts ... one has 
% the actual graphic ... the other (the bottom part) is a text
% table of all the attributes of interest for the graphic.

h0 = subplot('position',[0.1 0.3 0.8 0.6]); set(h0,'FontSize',14)
h = imagesc(plotdat.lons,plotdat.lats,datmat);
set(gca,'YDir','normal');
axis image;
worldmap;
h1 = colorbar;
set(h1,'FontSize',14)
set(get(h1,'Ylabel'),'String',plotdat.varunits,'FontSize',14)
title(sprintf('%s ''%s''',plotdat.model,var))

t0 = subplot('position',[0.1 0.05 0.8 0.25]);
t1 = plot([0 1],[0 1],'.','Visible','off');
axis off
text(0.0, 0.6, sprintf('%s %s',plotdat.model,fname),'Interpreter','none');
text(0.0, 0.4, timestring);
text(0.0, 0.2, sprintf('%s %s %s', copystring, var, plotdat.varlname));
text(0.0, 0.0, sprintf('level index %d at timestep %d', level, time));


% helper function(s) follow ...

function plotdat = netCDFinfo(fname,var)

ft                 = netcdf(fname);
plotdat.model      = ft.model(:);
plotdat.varlname   = ft{var}.long_name(:);
plotdat.varunits   = ft{var}.units(:);
plotdat.timeunits  = ft{'time'}.units(:);
plotdat.calendar   = ft{'time'}.calendar(:);
plotdat.dimcells   = dim(ft{var});
plotdat.dimnames   = ncnames(plotdat.dimcells);

for i = 1:length(plotdat.dimcells)

   dimstring = lower(ft{plotdat.dimnames{i}}.units(:));

   switch dimstring

      case 'degrees_east'
         plotdat.londimname = plotdat.dimnames{i};
         plotdat.lons       = getnc(fname,plotdat.dimnames{i});
         plotdat.lonunits   = dimstring;
      case 'degrees_north' 
         plotdat.latdimname = plotdat.dimnames{i};
         plotdat.lats       = getnc(fname,plotdat.dimnames{i});
         plotdat.latunits   = dimstring;
      case 'levels'
         plotdat.levdimname = plotdat.dimnames{i};
         plotdat.levels     = getnc(fname,plotdat.dimnames{i});
         plotdat.levunits   = dimstring;

      otherwise

   end

end
close(ft)

% Surface variables generally do not have a 'level'

if ( ~exist('plotdat.levdimname','var') )
   plotdat.levdimname = 'no level dimension';
   plotdat.levels     = 1;
   plotdat.levunits   = 'no level units';
end

