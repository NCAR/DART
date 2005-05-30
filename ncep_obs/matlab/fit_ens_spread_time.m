function fit_ens_spread_time(ddir)
% fit_ens_spread_time(ddir)
%
%
% USAGE
%
% fit_ens_spread_time('plot')

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

   varnames = {'T','W','Q','P'};

   Regions = {'Northern Hemisphere', ...
              'Southern Hemisphere', ...
              'Tropics', 'North America'};
   ptypes = {'gs-','bd-','ro-','k+-'};    % for each region

else
   error(sprintf('%s cannot be found.', datafile))
end

% set up a structure with all static plotting components

plotdat.toff      = toff;
plotdat.level     = level;
plotdat.pmax      = pmax;
plotdat.pmin      = pmin;
plotdat.linewidth = 2.0;
plotdat.flavor    = 'Ens Spread';
plotdat.ylabel    = 'RMSE'; 

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:length(varnames),

   % Set up a structure with all the plotting components

   plotdat.varname = varnames{ivar};

   switch obs_select
      case 1,
         obsstring = sprintf('%s (all data)',     plotdat.varname);
      case 2,
         obsstring = sprintf('%s (RaObs)',        plotdat.varname);
      otherwise,
         obsstring = sprintf('%s (ACARS,SATWND)', plotdat.varname);
   end

   switch varnames{ivar}
      case{'P'}
         ges  = sprintf('%sges_times.dat',varnames{ivar});
         anl  = sprintf('%sanl_times.dat',varnames{ivar});
         main = sprintf('%s %s',plotdat.flavor,obsstring);
      otherwise
         ges  = sprintf('%sges_times_%04dmb.dat',varnames{ivar},level);
         anl  = sprintf('%sanl_times_%04dmb.dat',varnames{ivar},level);
         main = sprintf('%s %s %d hPa',plotdat.flavor,obsstring,plotdat.level);
   end

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
%
subplot(2,2,plotdat.region)
plot(xp, yp_spread, 'k+-', xa, ya_spread, 'ro-', 'LineWidth', 1.5)
grid
ylabel(plotdat.ylabel, 'fontsize', 10)
ax = axis; ax(3) = 0.0;
axis(ax)
datetick('x',1)
title(plotdat.title, 'fontsize', 12, 'FontWeight', 'bold')
legend('guess', 'analysis')
h = legend;
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
