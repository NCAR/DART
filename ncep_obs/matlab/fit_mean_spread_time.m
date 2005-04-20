function fit_mean_spread_time(ddir)
% fit_mean_spread_time(ddir)
%
% ddir     is an optional argument specifying the directory containing
%               the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
%
% ddir = 'plot';
% fit_mean_spread_time(ddir)

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

% This ensures the directory with the datafiles 
% is in Matlab's search path.

if ( nargin > 0 )
   datafile = fullfile(ddir,'Tanl_times_level.dat');
else
   datafile = 'Tanl_times_level.dat';
   ddir = [];
end	

if ( exist(datafile,'file') == 2 )
   p     = load(datafile);
   level = p(1);
else
   error(sprintf('%s cannot be found.', datafile))
end

%----------------------------------------------------------------------
figure(1); clf; % Temperature
%----------------------------------------------------------------------

   % Set up a structure with all the plotting components
   plotdat.Regions = {'Northern Hemisphere', ...
                      'Southern Hemisphere', ...
                      'Tropics', 'North America'};
   plotdat.level   = level;
   plotdat.flavor  = 'ens mean';
   plotdat.xlabel  = 'Time interval';
   plotdat.ylabel  = 'RMSE';
   plotdat.ges     = fullfile(ddir,sprintf('Tges_times_%04dmb.dat',level));
   plotdat.anl     = fullfile(ddir,sprintf('Tanl_times_%04dmb.dat',level));
   plotdat.varname = 'T';

   plotdat.region = 1; myplot(plotdat);
   plotdat.region = 2; myplot(plotdat);
   plotdat.region = 3; myplot(plotdat);
   plotdat.region = 4; myplot(plotdat);

%----------------------------------------------------------------------
figure(2); clf; % Winds
%----------------------------------------------------------------------

   plotdat.ges     = fullfile(ddir,sprintf('Wges_times_%04dmb.dat',level));
   plotdat.anl     = fullfile(ddir,sprintf('Wanl_times_%04dmb.dat',level));
   plotdat.varname = 'Wind';

   plotdat.region = 1; myplot(plotdat);
   plotdat.region = 2; myplot(plotdat);
   plotdat.region = 3; myplot(plotdat);
   plotdat.region = 4; myplot(plotdat);

%----------------------------------------------------------------------
% common
%----------------------------------------------------------------------

print(1,'-dpsc','t_mean_spread_time.ps');
print(2,'-dpsc','w_mean_spread_time.ps');



function myplot(plotdat)
p  = load(plotdat.ges);
a  = load(plotdat.anl);

x = [1:2*size(p,1)];
ens_mean = x;
ens_spread = x;

countm = 2+(plotdat.region-1)*3;
counts = 3+(plotdat.region-1)*3;

for itime = 1:size(p,1)
   x(2*itime-1) = p(itime,1);
   ens_mean(2*itime-1) = p(itime,countm);
   ens_spread(2*itime-1) = p(itime,counts);
   x(2*itime) = p(itime,1);
   ens_mean(2*itime) = a(itime,countm);
   ens_spread(2*itime) = a(itime,counts);
end

subplot(2,2,plotdat.region)
   plot(x,ens_mean,'k+-',x,ens_spread,'ro-','LineWidth',1.5)
   grid
   xlabel(plotdat.xlabel, 'fontsize', 10);
   ylabel(plotdat.ylabel, 'fontsize', 10);
   string0 = sprintf('%s',plotdat.Regions{plotdat.region});
   string1 = sprintf('%s fit to RAobs', plotdat.varname);
   string2 = sprintf('%s %d hPa',plotdat.flavor,plotdat.level);
   title({string0,string1,string2}, 'fontsize', 12,'FontWeight','bold')
   h = legend('Ens. mean', 'Ens. spread');
   legend(h,'boxoff')
