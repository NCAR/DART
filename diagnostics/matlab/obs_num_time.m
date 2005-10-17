function obs_num_time(ddir)
% obs_num_time(ddir)
%
% Part of the observation-space diagnostics routines.
%
% Plots the total number of observations as a function of time 
% for a given level for several different regions.
% 
% ddir   an optional argument specifying the directory containing
%        the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
% 
% obs_num_time('plot')
%
% USAGE: if the preprocessed data files are in the current directory
%
% obs_num_time


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

datafile = 'ObsDiagAtts';
ptypes = {'gs-','bd-','ro-','k+-'};    % for each region

%----------------------------------------------------------------------
% Get attributes from obs_diag run.
%----------------------------------------------------------------------

if ( exist(datafile) == 2 )

   eval(datafile)

   if ( exist('plevel','var') == 0 )
      plevel = 1;
   end

else
   error(sprintf('%s cannot be found.', datafile))
end

% set up a structure with all static plotting components
temp = datenum(obs_year,obs_month,obs_day);
plotdat.toff      = temp - round(t1); % determine temporal offset (calendar base)
plotdat.day1      = datestr(t1+plotdat.toff,'yyyy-mm-dd HH');
plotdat.dayN      = datestr(tN+plotdat.toff,'yyyy-mm-dd HH');
plotdat.obs_year  = obs_year;
plotdat.obs_month = obs_month;
plotdat.obs_day   = obs_day;
plotdat.level     = plevel;
plotdat.ylabel    = 'observation count';
plotdat.nregions  = length(Regions);
plotdat.nvars     = length(One_Level_Varnames);

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars,

   % set up a structure with all the plotting components

   plotdat.varname = One_Level_Varnames{ivar};

   switch obs_select
      case 1,
         string1 = sprintf('%s (all data)',     plotdat.varname);
      case 2, 
         string1 = sprintf('%s (RaObs)',        plotdat.varname);
      otherwise,
         string1 = sprintf('%s (ACARS,SATWND)', plotdat.varname);
   end

%  switch One_Level_Varnames{ivar}
%     case{'P'}
         plotdat.fname = sprintf('%s_ges_times.dat',One_Level_Varnames{ivar});
         plotdat.main  = sprintf('%s',string1);
%     otherwise
%        plotdat.fname = sprintf('%s_ges_times_%04dmb.dat',One_Level_Varnames{ivar},plotdat.level);
%        plotdat.main  = sprintf('%s %d hPa',string1,plotdat.level);
%  end

   % plot by region

   page1 = 2*(ivar-1)+1;
   page2 = 2*(ivar-1)+2;
   figure(page1); clf;

   for iregion = 1:plotdat.nregions,
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

   if (plotdat.nregions > 3 )
      h = plot(xax, NbyRegion(:,1), ptypes{1}, ...
               xax, NbyRegion(:,2), ptypes{2}, ...
	       xax, NbyRegion(:,3), ptypes{3}, ...
	       xax, NbyRegion(:,4), ptypes{4}, 'LineWidth', 2.0);
      h = legend(Regions{1},Regions{2},Regions{3},Regions{4}, ...
	         'Location','NorthWest');
   else
      h = plot(xax, NbyRegion(:,1), ptypes{1}, ...
               xax, NbyRegion(:,2), ptypes{2}, ...
	       xax, NbyRegion(:,3), ptypes{3}, 'LineWidth', 2.0);
      h = legend(Regions{1},Regions{2},Regions{3}, ...
   	         'Location','NorthWest');
   end

   grid
   title(plotdat.main,'Interpreter', 'none', ...
       'FontSize', 12, 'FontWeight', 'bold')
   ylabel(plotdat.ylabel, 'fontsize', 10)
   datetick('x',1)
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
   xlabel('days')

   if (plotdat.obs_year > 1000); datetick('x',1); end
   
   datetick('x',1)
   ylabel(plotdat.ylabel, 'FontSize', 10) ;
   title(plotdat.title,'Interpreter', 'none', ...
        'FontSize', 12, 'FontWeight', 'bold')
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
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','bottom', ...
      'Interpreter','none', ...
      'FontSize',12, ...
      'FontWeight','bold')



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
