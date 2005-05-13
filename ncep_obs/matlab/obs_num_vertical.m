function obs_num_vertical(ddir)
% obs_num_vertical    Plots the number of observations as a function of height for several different regions.
% 
% ddir     is an optional argument specifying the directory containing
%          the data files as preprocessed by the support routines.
% 
% USAGE: if the preprocessed data files are in a directory called 'plot'
% 
% obs_num_vertical('plot')

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
   day1 = datestr(t1+toff+iskip,'yyyy-mm-dd HH');
   dayN = datestr(tN+toff,      'yyyy-mm-dd HH');
   pmax = psurface;
   pmin = ptop;

   % There is no vertical distribution of surface pressure

   varnames = {'T','W','Q'};

   Regions = {'Northern Hemisphere', ...
              'Southern Hemisphere', ...
              'Tropics', 'North America'};
   ptypes = {'gs-','bd-','ro-','k+-'};    % for each region

else
   error(sprintf('%s cannot be found.', datafile))
end

% set up a structure with all static plotting components

plotdat.toff      = toff;
plotdat.linewidth = 2.0;
plotdat.pmax      = pmax;
plotdat.pmin      = pmin;
plotdat.ylabel    = 'Pressure (hPa)';
plotdat.xlabel    = 'observation count';

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

   plotdat.fname = sprintf('%sges_ver_ave.dat',varnames{ivar});
   plotdat.main  = sprintf('%s %s -- %s',string1,day1,dayN);

   % plot by region

   page1 = 2*(ivar-1)+1;
   page2 = 2*(ivar-1)+2;
   figure(page1); clf;

   for iregion = 1:length(Regions),
      plotdat.title  = Regions{iregion};
      plotdat.region = iregion;
      plotdat.ptype  = ptypes{iregion};
      [nobs,obslevels] = myplot(plotdat);
      NbyRegion(:,iregion) = nobs;
   end

   CenterAnnotation(plotdat.main)
   BottomAnnotation(plotdat.fname)

   psfname = sprintf('%s_num_vert.ps',plotdat.varname);
   print(page1,'-dpsc',psfname);

   % All regions on one figure

   figure(page2); clf;

   h = plot(NbyRegion(:,1), obslevels, ptypes{1}, ...
            NbyRegion(:,2), obslevels, ptypes{2}, ...
            NbyRegion(:,3), obslevels, ptypes{3}, ...
            NbyRegion(:,4), obslevels, ptypes{4}, 'LineWidth', plotdat.linewidth);
   grid
   ax = axis; 
   ax(3) = plotdat.pmin; 
   ax(4) = plotdat.pmax; 
   axis(ax)
   set(gca,'YDir', 'reverse')
   title(plotdat.main, 'FontSize', 12, 'FontWeight', 'bold')
   ylabel(plotdat.ylabel, 'fontsize', 10)
   xlabel(plotdat.xlabel, 'fontsize', 10)
   
   h = legend(Regions{1},Regions{2},Regions{3},Regions{4}, ...
              'Location','Best');
   legend(h,'boxoff');

   BottomAnnotation(plotdat.fname)

   str = sprintf('print -f%d -dpsc -append %s',page2,psfname);
   eval(str)

end

path(startpath); % restore MATLABPATH to original setting

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------
function [Nobs,levels] = myplot(plotdat)

if ( exist(plotdat.fname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',plotdat.fname))
end

regionindex = 3 + 2*(plotdat.region - 1);

datmat = load(plotdat.fname);
obsmat = SqueezeMissing(datmat);
levels = obsmat(:,1); 
Nobs   = obsmat(:,regionindex); 

subplot(2,2,plotdat.region)
   plot(Nobs, levels, plotdat.ptype , 'LineWidth', plotdat.linewidth)
   grid
   ax = axis; 
   ax(3) = plotdat.pmin; 
   ax(4) = plotdat.pmax; 
   axis(ax)
   set(gca,'YDir', 'reverse')
   title( plotdat.title,  'fontsize', 12,'FontWeight','bold')
   ylabel(plotdat.ylabel, 'fontsize', 10)
   xlabel(plotdat.xlabel, 'fontsize', 10)



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
      'FontSize',8)
