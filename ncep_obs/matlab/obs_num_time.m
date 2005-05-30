function obs_num_time(ddir)
% obs_num_time     Plots the total number of observations as a function of time for a given level for several different regions.
% 
% ddir     an optional argument specifying the directory containing
%          the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
% 
% obs_num_time('plot')

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

datafile = 'ObsDiagAtts';

%----------------------------------------------------------------------
% Get attributes from obs_diag run.
%----------------------------------------------------------------------

if ( exist(datafile) == 2 )

   eval(datafile)

   temp = datenum(obs_year,obs_month,obs_day);
   toff = temp - round(t1); % determine temporal offset (calendar base)
   day1 = datestr(t1+toff,'yyyy-mm-dd HH');
   dayN = datestr(tN+toff,'yyyy-mm-dd HH');
   pmax = psurface;
   pmin = ptop;

   % There is no vertical distribution of surface pressure

   varnames = {'T','W','Q'};
   varnames = {'T','W','Q','P'};

   Regions = {'Northern Hemisphere', ...
              'Southern Hemisphere', ...
              'Tropics', 'North America'};
   ptypes = {'gs-','bd-','ro-','k+-'};    % for each region

else
   error(sprintf('%s cannot be found.', datafile))
end

% set up a structure with all static plotting components

plotdat.level   = level;
plotdat.toff    = toff;
plotdat.ylabel  = 'observation count';

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:length(varnames),

   % set up a structure with all the plotting components

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
         plotdat.fname = sprintf('%sges_times.dat',varnames{ivar});
         plotdat.main  = sprintf('%s',string1);
      otherwise
         plotdat.fname = sprintf('%sges_times_%04dmb.dat',varnames{ivar},level);
         plotdat.main  = sprintf('%s %d hPa',string1,plotdat.level);
   end

   % plot by region

   page1 = 2*(ivar-1)+1;
   page2 = 2*(ivar-1)+2;
   figure(page1); clf;

   for iregion = 1:length(Regions),
      plotdat.title  = Regions{iregion};
      plotdat.region = iregion;
      plotdat.ptype  = ptypes{iregion};
      [nT,xax]       = myplot(plotdat);
      NbyRegion(:,iregion) = nT;
   end

   CenterAnnotation(plotdat.main)
   BottomAnnotation(plotdat.fname)

   psfname = sprintf('%s_obs_num_time.ps',plotdat.varname);
   print(page1,'-dpsc',psfname);

   % All regions on one figure

   figure(page2); clf;

   h = plot(xax, NbyRegion(:,1), ptypes{1}, ...
            xax, NbyRegion(:,2), ptypes{2}, ...
	    xax, NbyRegion(:,3), ptypes{3}, ...
	    xax, NbyRegion(:,4), ptypes{4}, 'LineWidth', 2.0);
   grid
   title(plotdat.main, 'FontSize', 12, 'FontWeight', 'bold')
   ylabel(plotdat.ylabel, 'fontsize', 10)
   datetick('x',1)
   h = legend(Regions{1},Regions{2},Regions{3},Regions{4}, ...
	      'Location','NorthWest');
   legend(h,'boxoff');

   BottomAnnotation(plotdat.fname)

   str = sprintf('print -f%d -dpsc -append %s',page2,psfname);
   eval(str)

end

path(startpath); % restore MATLABPATH to original setting

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------

function [yp_num,xp] = myplot(plotdat)

p1 = load(plotdat.fname); p = SqueezeMissing(p1);
xp = p(:,1) + p(:,2)/86400 + plotdat.toff;

offset = 5;  % columns 1,2 are time, 3=mean, 4=spread, 5=numobs

count  = offset+(plotdat.region-1)*3;
yp_num = p(:,count);

subplot(2,2,plotdat.region)
   plot(xp,yp_num,plotdat.ptype,'LineWidth',2.0)
   grid
   datetick('x',1)
   ylabel(plotdat.ylabel, 'fontsize', 10) ;
   title(plotdat.title, 'fontsize', 12,'FontWeight','bold')
   axis([min(xp) max(xp) -Inf Inf])



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
