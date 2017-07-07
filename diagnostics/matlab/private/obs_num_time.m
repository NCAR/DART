function obs_num_time(ddir)
%% obs_num_time(ddir)
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

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

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

if ( exist(datafile,'file') == 2 )

   eval(datafile)

   if ( exist('plevel','var') == 0 )
      plevel = 1;
      iskip = iskip_days;
      plotdat.toff = 0;
      plotdat.bin1 = datenum(t1);
   else  % high dimensional models
      % Coordinate between time types and dates
      skip_seconds = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
      iskip = time_to_skip(3) + skip_seconds/86400;
      plotdat.bin1 = datenum(first_bin_center);
      plotdat.toff = plotdat.bin1 - t1;
   end

else
   error('%s cannot be found.', datafile)
end

% set up a structure with all static plotting components
plotdat.day1      = datestr(t1+plotdat.toff+iskip,'yyyy-mm-dd HH');
plotdat.dayN      = datestr(tN+plotdat.toff,'yyyy-mm-dd HH');
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
   figure(page1); orient tall; clf; wysiwyg

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
   % We assume there is at least one region and then we'll add to that.

   figure(page2); clf;

   h = plot( xax, NbyRegion(:,1), ptypes{1}, 'LineWidth', 2.0);
   [~, ~, outh, outm] = legend(Regions{1});
   hold on;

   for iregion = 2:length(Regions),
      h = plot( xax, NbyRegion(:,iregion), ptypes{iregion}, 'LineWidth', 2.0);
      % Must add salient information to the legend.
      % legh     handle to the legend axes
      % objh     handle for the text, lines, and patches in the legend
      % outh     handle for the lines and patches in the plot
      % outm     cell array for the text in the legend
      nlines = length(outm);
      outm{nlines + 1} = Regions{iregion};
      [~, ~, outh, outm] = legend([outh; h],outm,'Location','NorthEast');
   end

   legend boxoff
   grid
   title(plotdat.main,'Interpreter', 'none', ...
       'FontSize', 12, 'FontWeight', 'bold')
   ylabel(plotdat.ylabel, 'fontsize', 10)

   % a slightly better way to annotate dates, etc.
   ttot = max(xax) - min(xax) + 1;
   if ((plotdat.bin1 > 1000) && (ttot > 32));
      datetick('x',6,'keeplimits','keepticks');
      monstr = datestr(xax(1),28);
      xlabel(sprintf('month/day - %s start',monstr))
   elseif (plotdat.bin1 > 1000);
      datetick('x',7);
      monstr = datestr(xax(1),28);
      xlabel(sprintf('day of month - %s start',monstr))
   else
      xlabel('days')
   end

   BottomAnnotation(plotdat.fname)

   str = sprintf('print -f%d -dpsc -append %s',page2,psfname);
   eval(str)
   hold off;

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


subplot(plotdat.nregions,1,plotdat.region)

   plot(xp,yp_num,plotdat.ptype,'LineWidth',2.0)
   axis([min(xp) max(xp) -Inf Inf])
   grid
   ylabel(plotdat.ylabel, 'FontSize', 10) ;
   title(plotdat.title,'Interpreter', 'none', ...
        'FontSize', 12, 'FontWeight', 'bold')

   % a slightly better way to annotate dates, etc.
   ttot = max(xp) - min(xp) + 1;
   if ((plotdat.bin1 > 1000) && (ttot > 32));
      datetick('x',6,'keeplimits','keepticks');
      monstr = datestr(xp(1),28);
      xlabel(sprintf('month/day - %s start',monstr))
   elseif (plotdat.bin1 > 1000);
      datetick('x',7);
      monstr = datestr(xp(1),28);
      xlabel(sprintf('day of month - %s start',monstr))
   else
      xlabel('days')
   end



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
[pathstr,~,~] = fileparts(bob);
h = text(0.0,0.5,pathstr);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',8)


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
