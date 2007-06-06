function fit_ens_spread_time(ddir)
% fit_ens_spread_time(ddir)
%
% Part of the observation-space diagnostics routines.
%
% Plots the spatial mean RMS of the spread of the ensemble as a function 
% of time for both the 'guess' and the 'analysis' at a single level. 
% Several regions are plotted. This function simply plots the 
% data in *ges_times.dat using metadata in ObsDiagAtts.m - both
% created by the executable 'obs_diag'.
%
% 'obs_diag' also produces a matlab-compatible file of plotting attributes:
% ObsDiagAtts.m which specifies the run-time configuration of obs_diag.
%
% The figures are automatically saved as postscript files.
%
% ddir   is an optional argument specifying the directory containing
%        the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
%
% fit_ens_spread_time('plot')
%
% USAGE: if the preprocessed data files are in the current directory
%
% fit_ens_spread_time

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
      % set up a structure with all static plotting components
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
plotdat.ylabel    = 'RMS'; 
plotdat.nregions  = length(Regions); 
plotdat.nvars     = length(One_Level_Varnames); 
plotdat.flavor    = 'Ens Spread';
plotdat.linewidth = 2.0;

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars,

   plotdat.varname = One_Level_Varnames{ivar};

   switch obs_select
      case 1,
         obsstring = sprintf('%s (all data)',     plotdat.varname);
      case 2,
         obsstring = sprintf('%s (RaObs)',        plotdat.varname);
      otherwise,
         obsstring = sprintf('%s (ACARS,SATWND)', plotdat.varname);
   end

%  switch One_Level_Varnames{ivar}
%     case{'P'}
         ges  = sprintf('%s_ges_times.dat',One_Level_Varnames{ivar});
         anl  = sprintf('%s_anl_times.dat',One_Level_Varnames{ivar});
         main = sprintf('%s %s',plotdat.flavor,obsstring);
%     otherwise
%        ges  = sprintf('%s_ges_times_%04dmb.dat',One_Level_Varnames{ivar},plotdat.level);
%        anl  = sprintf('%s_anl_times_%04dmb.dat',One_Level_Varnames{ivar},plotdat.level);
%        main = sprintf('%s %s %d hPa',plotdat.flavor,obsstring,plotdat.level);
%  end

   plotdat.ges     = ges;
   plotdat.anl     = anl;

   % plot each region

   figure(ivar); clf;

   for iregion = 1:length(Regions),
      plotdat.title   = Regions{iregion};
      plotdat.region  = iregion;
      Myplot(plotdat)
   end

   CenterAnnotation(main);  % One title in the middle
   BottomAnnotation(ges);   % directory in middle, bottom

   % create a postscript file

   psfname = sprintf('%s_ens_spread_time.ps',plotdat.varname);
   print(ivar,'-dpsc',psfname);

end

path(startpath); % restore MATLABPATH to original setting

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------


function Myplot(plotdat)
%
% After the first column, each set of 3 columns
% represents a different region.
% Depends on the format written by obs_diag.f90
%
p1   = load(plotdat.ges); p = SqueezeMissing(p1);
a1   = load(plotdat.anl); a = SqueezeMissing(a1);

xp        = p(:,1) + p(:,2)./86400 + plotdat.toff;
xa        = a(:,1) + a(:,2)./86400 + plotdat.toff;

offset    = 4;  % columns 1,2 are time, 3=mean, 4=spread, 5=numobs

count     = offset+(plotdat.region-1)*3;
yp_spread = p(:,count);
ya_spread = a(:,count);

gmean = mean(yp_spread(isfinite(yp_spread))); gstring = sprintf('guess;    mean=%.3f',gmean);
amean = mean(ya_spread(isfinite(ya_spread))); astring = sprintf('analysis; mean=%.3f',amean);

subplot(2,2,plotdat.region)
   plot(xp, yp_spread, 'k+-', xa, ya_spread, 'ro-', 'LineWidth', 1.5)
   grid
   ax = axis; ax(3) = 0.0; axis(ax)

   if (plotdat.bin1 > 1000); 
      datetick('x',1);
   else
      xlabel('days')
   end

   ylabel(plotdat.ylabel, 'fontsize', 10)
   title(plotdat.title,'Interpreter', 'none', 'fontsize', 12, 'FontWeight', 'bold')
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
