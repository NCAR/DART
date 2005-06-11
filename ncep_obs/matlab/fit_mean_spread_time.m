function fit_mean_spread_time(ddir)
% fit_mean_spread_time(ddir)
%
% Plots the Ensemble mean and spread as a function of time at a single 
% level for several regions. This function simply plots the 
% data in *ges_times_*mb.dat using metadata in ObsDiagAtts.m - both
% created by the executable 'obs_diag'.

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

% Ensures the specified directory is searched first.
if ( nargin > 0 )
   startpath = addpath(ddir);
else
   startpath = path;
end

%----------------------------------------------------------------------
% Defaults
%----------------------------------------------------------------------

varnames = {'T','W','Q','P'};

Regions = {'Northern Hemisphere', ...
           'Southern Hemisphere', ...
           'Tropics', 'North America'};

ptypes = {'gs-','bd-','ro-','k+-'};    % for each region

%----------------------------------------------------------------------
% Get attributes from obs_diag run.
%----------------------------------------------------------------------

datafile = 'ObsDiagAtts';

if ( exist(datafile) == 2 )

   eval(datafile)

   temp = datenum(obs_year,obs_month,obs_day);
   toff = temp - round(t1); % determine temporal offset (calendar base)
   day1 = datestr(t1+toff,'yyyy-mm-dd HH');
   dayN = datestr(tN+toff,'yyyy-mm-dd HH');

else
   error(sprintf('%s cannot be found.', datafile))
end

% set up a structure with all static plotting components

plotdat.level     = level;
plotdat.toff      = toff;
plotdat.ylabel    = 'RMSE';
plotdat.nregions  = length(Regions);
plotdat.nvars     = length(varnames);
plotdat.linewidth = 2.0;

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars,

   plotdat.varname = varnames{ivar};

   switch obs_select
      case 1,
         string1 = sprintf('%s (all data)',     plotdat.varname);
      case 2, 
         string1 = sprintf('%s (RaObs)',        plotdat.varname);
      otherwise,
         string1 = sprintf('%s (ACARS,SATWND)', plotdat.varname);
   end

   switch varnames{ivar}
      case{'P'}
         ges  = sprintf('%sges_times.dat',varnames{ivar});
         anl  = sprintf('%sanl_times.dat',varnames{ivar});
         main = sprintf('%s',string1);
      otherwise
         ges  = sprintf('%sges_times_%04dmb.dat',varnames{ivar},level);
         anl  = sprintf('%sanl_times_%04dmb.dat',varnames{ivar},level);
         main = sprintf('%s %d hPa',string1,plotdat.level);
   end

   plotdat.ges     = ges;
   plotdat.anl     = anl;

   % plot by region

   figure(ivar); clf;

   for iregion = 1:plotdat.nregions,
      plotdat.title  = Regions{iregion};
      plotdat.region = iregion;
      myplot(plotdat);
   end

   CenterAnnotation(main)
   BottomAnnotation(ges)

   % create a postscript file

   psfname = sprintf('%s_mean_spread_time.ps',plotdat.varname);
   print(ivar,'-dpsc',psfname);

end

path(startpath); % restore MATLABPATH to original setting

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------

function myplot(plotdat)
p1 = load(plotdat.ges); p = SqueezeMissing(p1);
a1 = load(plotdat.anl); a = SqueezeMissing(a1);

x = [1:2*size(p,1)]; % Each time has a guess and analysis 
ens_mean   = x;
ens_spread = x;

countm = 3+(plotdat.region-1)*3; % pick off region mean
counts = 4+(plotdat.region-1)*3; % pick off region spread

for itime = 1:size(p,1)

   obsT1 = p(itime,1) + p(itime,2)/86400 + plotdat.toff;
   obsT2 = a(itime,1) + a(itime,2)/86400 + plotdat.toff;
   x(2*itime-1) = obsT1;
   x(2*itime  ) = obsT2;

   ens_mean(  2*itime-1) = p(itime,countm);
   ens_spread(2*itime-1) = p(itime,counts);

   ens_mean(  2*itime  ) = a(itime,countm);
   ens_spread(2*itime  ) = a(itime,counts);
end

subplot(2,2,plotdat.region)
   % Since the mean and the spread are getting plotted
   % on the same figure, we should have two axes ... 
   % bias on left, spread on right, for example. no time now ...

   gmean = mean(ens_mean);   gstring = sprintf('Ens. mean;   mean=%.3f',gmean);
   amean = mean(ens_spread); astring = sprintf('Ens. spread; mean=%.3f',amean);

   plot(x,ens_mean,'k+-',x,ens_spread,'ro-','LineWidth',1.5)
   grid
   ylabel(plotdat.ylabel, 'fontsize', 10);
   datetick('x',1)
   title(plotdat.title, 'fontsize', 12,'FontWeight','bold')
   h = legend(gstring, astring);
   legend(h,'boxoff')


function y = SqueezeMissing(x)

missing = find(x < -98); % 'missing' is coded as -99

if isempty(missing)
  y = x;
else
  y = x;
  y(missing) = NaN;
end



function CenterAnnotation(main)
subplot('position',[0.48 0.48 0.04 0.04])
axis off
h = text(0.5,0.5,main);
set(h,'HorizontalAlignment','center','VerticalAlignment','bottom',...
   'FontSize',12,'FontWeight','bold')



function BottomAnnotation(main)
% annotates the directory containing the data being plotted
subplot('position',[0.48 0.01 0.04 0.04])
axis off
bob = which(main);
[pathstr,name,ext,versn] = fileparts(bob);
h = text(0.0,0.5,pathstr);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',8)
