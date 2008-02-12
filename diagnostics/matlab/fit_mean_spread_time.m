function fit_mean_spread_time(ddir)
% fit_mean_spread_time(ddir)
%
% Part of the observation-space diagnostics routines.
%
% Plots the spatial mean RMSE of the ensemble mean and 
% RMS of the ensemble spread as a function of time at a single 
% level for several regions. This function simply plots the 
% data in *ges_times_*mb.dat using metadata in ObsDiagAtts.m - both
% created by the executable 'obs_diag'.
%
% 'obs_diag' also produces a matlab-compatible file of plotting attributes:
% ObsDiagAtts.m which specifies the run-time configuration of obs_diag.
%
% ddir     is an optional argument specifying the directory containing
%               the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
%
% ddir = 'plot';
% fit_mean_spread_time(ddir)
%
% USAGE: if the preprocessed data files are in the current directory
%
% fit_mean_spread_time

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
% Get plotting metadata from obs_diag run.
%----------------------------------------------------------------------

if ( exist(datafile) == 2 )

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
      plotdat.bin1 = datenum(first_bin_center); % a known date in matlab's time units
      plotdat.toff = plotdat.bin1 - t1;         % determine temporal offset (calendar base)
   end

else
   error(sprintf('%s cannot be found.', datafile))
end

% set up a structure with all static plotting components
plotdat.day1      = datestr(t1+plotdat.toff+iskip,'yyyy-mm-dd HH');
plotdat.dayN      = datestr(tN+plotdat.toff,'yyyy-mm-dd HH');
plotdat.level     = plevel; 
plotdat.ylabel    = 'RMSE';
plotdat.nregions  = length(Regions);
plotdat.nvars     = length(One_Level_Varnames);
plotdat.linewidth = 2.0;

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars,

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
         ges  = sprintf('%s_ges_times.dat',One_Level_Varnames{ivar});
         anl  = sprintf('%s_anl_times.dat',One_Level_Varnames{ivar});
         main = sprintf('%s',string1);
%     otherwise
%        ges  = sprintf('%s_ges_times_%04dmb.dat',One_Level_Varnames{ivar},plotdat.level);
%        anl  = sprintf('%s_anl_times_%04dmb.dat',One_Level_Varnames{ivar},plotdat.level);
%        main = sprintf('%s %d hPa',string1,plotdat.level);
%  end

   plotdat.ges     = ges;
   plotdat.anl     = anl;

   % plot by region

   figure(ivar); orient tall; clf; wysiwyg

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

subplot(plotdat.nregions,1,plotdat.region)

   % Since the mean and the spread are getting plotted
   % on the same figure, we should have two axes ... 
   % bias on left, spread on right, for example. no time now ...

   gmean = mean(ens_mean(  isfinite(ens_mean  )));   
   amean = mean(ens_spread(isfinite(ens_spread))); 
   gstring = sprintf('Ens. mean;   mean=%.3f',gmean);
   astring = sprintf('Ens. spread; mean=%.3f',amean);

   plot(x,ens_mean,'k+-',x,ens_spread,'ro-','LineWidth',1.5)
   grid
   ylabel(plotdat.ylabel, 'fontsize', 10);
   title(plotdat.title, 'Interpreter', 'none', ...
         'Fontsize', 12, 'FontWeight', 'bold')
   h = legend(gstring, astring);
   legend(h,'boxoff')

   % a slightly better way to annotate dates, etc.
   ttot = max(x) - min(x) + 1;
   if ((plotdat.bin1 > 1000) && (ttot > 32));
      datetick('x',6,'keeplimits','keepticks');
      monstr = datestr(x(1),28);
      xlabel(sprintf('month/day - %s start',monstr))
   elseif (plotdat.bin1 > 1000);
      datetick('x',7);
      monstr = datestr(x(1),28);
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
[pathstr,name,ext,versn] = fileparts(bob);
h = text(0.0,0.5,pathstr);
set(h,'HorizontalAlignment','center', ...
      'VerticalAlignment','middle',...
      'Interpreter','none',...
      'FontSize',8)
