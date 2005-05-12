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

% Make sure the directory with the files is available to matlab
if (nargin > 0 )
   datafile = fullfile(ddir,'Tanl_times_level');
else
   datafile = 'Tanl_times_level';
   ddir = [];
end

%----------------------------------------------------------------------
% Get attributes from obs_diag run.
%----------------------------------------------------------------------

if ( exist(datafile) == 2 )

   eval(datafile)

   temp   = datenum(obs_year,obs_month,obs_day);
   toff = temp - round(t1); % determine temporal offset (calendar base)
   day1 = datestr(t1+toff,'yyyy-mm-dd HH');
   dayN = datestr(tN+toff,'yyyy-mm-dd HH');

else
   error(sprintf('%s cannot be found.', datafile))
end

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

varnames = {'T','W','Q','P'};
Regions = {'Northern Hemisphere', ...
           'Southern Hemisphere', ...
           'Tropics', 'North America'};
ptypes = {'gs-','bd-','ro-','k+-'};    % for each region

for ivar = 1:length(varnames),

   % Set up a structure with all the plotting components

   plotdat.level   = level;
   plotdat.ylabel  = 'observation count';
   plotdat.varname = varnames{ivar};
   plotdat.toff    = toff;

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
         fname = fullfile(ddir,sprintf('%sges_times.dat',varnames{ivar}));
         main = sprintf('%s',string1);
      otherwise
         fname = fullfile(ddir,sprintf('%sges_times_%04dmb.dat',varnames{ivar},level));
         main = sprintf('%s %d hPa',string1,plotdat.level);
   end

   plotdat.fname = fname;

   % plot by region

   page1 = 2*(ivar-1)+1;
   page2 = 2*(ivar-1)+2;
   figure(page1); clf;

%  NbyRegion = zeros(xxx,iregion)
   for iregion = 1:length(Regions),
      plotdat.title  = Regions{iregion}; 
      plotdat.region = iregion;
      plotdat.ptype  = ptypes{iregion};
      [nT,xax]       = myplot(plotdat);
      NbyRegion(:,iregion) = nT;
   end

   %[nT_SH,xax]    = myplot(plotdat);
   %[nT_TR,xax]    = myplot(plotdat);
   %[nT_NA,xax]    = myplot(plotdat);

   CenterAnnotation(main)

   psfname = sprintf('%s_obs_num_time.ps',plotdat.varname);
   print(page1,'-dpsc',psfname);
  
   % All regions on one figure

   figure(page2); clf;

   h = plot(xax, NbyRegion(:,1), ptypes{1}, ...
            xax, NbyRegion(:,2), ptypes{2}, ...
	    xax, NbyRegion(:,3), ptypes{3}, ...
	    xax, NbyRegion(:,4), ptypes{4}, 'LineWidth', 2.0);
   grid
   title(main, 'FontSize', 12, 'FontWeight', 'bold')
   ylabel(plotdat.ylabel, 'fontsize', 10)
   datetick('x',1)
   h = legend(Regions{1},Regions{2},Regions{3},Regions{4}, ...
	      'Location','NorthWest');
   legend(h,'boxoff');

   str = sprintf('print -f%d -dpsc -append %s',page2,psfname);
   eval(str)

end

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
   title(plotdat.title, 'fontsize', 14,'FontWeight','bold')



function y = SqueezeMissing(x)

missing = find(x < -98); % 'missing' is coded as -99

if isempty(missing)
  y = x;
else
  y = x;
  y(missing) = NaN;
end


function CenterAnnotation(top)
subplot('position',[0.48 0.48 0.04 0.04])
axis off
h = text(0.5,0.5,top);
set(h,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
   'FontSize',12,'FontWeight','bold')
