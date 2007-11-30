function obs_num_vertical(ddir)
% obs_num_vertical(ddir)
%
% Plots the number of observations as a function of height for 
% several different regions.
% 
% ddir     is an optional argument specifying the directory containing
%          the data files as preprocessed by the support routines.
% 
% USAGE: if the preprocessed data files are in a directory called 'plot'
% 
% obs_num_vertical('plot')

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

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
      disp(sprintf('%s does not have multiple levels.', datafile))
      disp('It cannot be plotted with obs_num_vertical.')
      return
   end

else
   error(sprintf('%s cannot be found.', datafile))
end

% set up a structure with all static plotting components
skip_seconds = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
iskip = time_to_skip(3) + skip_seconds/86400;

plotdat.bin1      = datenum(first_bin_center); % a known date in matlab's time units
plotdat.toff      = plotdat.bin1 - t1;         % determine temporal offset (calendar base)
plotdat.day1      = datestr(t1+plotdat.toff+iskip,'yyyy-mm-dd HH');
plotdat.dayN      = datestr(tN+plotdat.toff,'yyyy-mm-dd HH');
plotdat.linewidth = 2.0;
plotdat.xlabel    = 'observation count';

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:length(All_Level_Varnames),

   % set up a structure with all the plotting components

   plotdat.varname = All_Level_Varnames{ivar};

   switch obs_select
      case 1,
         string1 = sprintf('%s (all data)',     plotdat.varname);
      case 2, 
         string1 = sprintf('%s (RaObs)',        plotdat.varname);
      otherwise,
         string1 = sprintf('%s (ACARS,SATWND)', plotdat.varname);
   end

   plotdat.fname = sprintf('%s_ges_ver_ave.dat',All_Level_Varnames{ivar});
   plotdat.main  = sprintf('%s %s -- %s',string1,plotdat.day1,plotdat.dayN);
   plotdat       = SetLevels(plotdat);

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
   % We assume there is at least one region and then we'll add to that.

   figure(page2); clf;

   h = plot(NbyRegion(:,1), obslevels, ptypes{1}, 'LineWidth', plotdat.linewidth);
   [legh, objh, outh, outm] = legend(Regions{1});
   hold on;

   for iregion = 2:length(Regions),
      h = plot(NbyRegion(:,iregion), obslevels, ptypes{iregion}, ...
                                'LineWidth', plotdat.linewidth);
      % Must add salient information to the legend.
      % legh     handle to the legend axes
      % objh     handle for the text, lines, and patches in the legend
      % outh     handle for the lines and patches in the plot
      % outm     cell array for the text in the legend
      nlines = length(outm);
      outm{nlines + 1} = Regions{iregion};
      [legh, objh, outh, outm] = legend([outh; h],outm,0);
      %h = legend(Regions{1},Regions{2},Regions{3},Regions{4}, 'Location','Best');
   end

   legend boxoff
   grid
   ax = axis; 
   axis([ax(1) ax(2) min(plotdat.ylims) max(plotdat.ylims)])
   set(gca,'YDir', plotdat.ydir)
   title(plotdat.main, 'Interpreter','none','FontSize', 12, 'FontWeight', 'bold')
   ylabel(plotdat.ylabel, 'fontsize', 10)
   xlabel(plotdat.xlabel, 'fontsize', 10)
   BottomAnnotation(plotdat.fname)

   str = sprintf('print -f%d -dpsc -append %s',page2,psfname);
   eval(str)

   hold off

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
levels = datmat(:,1); 
obsmat = SqueezeMissing(datmat);
Nobs   = obsmat(:,regionindex); 

subplot(2,2,plotdat.region)
   plot(Nobs, levels, plotdat.ptype , 'LineWidth', plotdat.linewidth)
   grid
   ax = axis; 
   axis([ax(1) ax(2) min(plotdat.ylims) max(plotdat.ylims)])
   set(gca,'YDir', plotdat.ydir)

   title( plotdat.title,  'Interpreter','none','fontsize', 12,'FontWeight','bold')
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



function plotstruct = SetLevels(plotdat)

if ( exist(plotdat.fname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',plotdat.fname))
end

plotstruct = plotdat;

datmat = load(plotstruct.fname);
levels = datmat(:,1); 

plotstruct.top     = levels(1);
plotstruct.surface = levels(length(levels));

if ( plotstruct.top > plotstruct.surface )
   plotstruct.ylabel = 'height(m)';
   plotstruct.ydir   = 'normal';
   plotstruct.ylims  = [plotstruct.surface plotstruct.top];
else
   plotstruct.ylabel = 'Pressure (hPa)';
   plotstruct.ydir   = 'reverse';
   plotstruct.ylims  = [plotstruct.surface plotstruct.top];
end
