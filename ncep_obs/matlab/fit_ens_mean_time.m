function fit_ens_mean_time(ddir)
% fit_ens_mean_time(ddir)
%
% fit_ens_mean_time: plots the average of all the observations in given regions
%          as a function of time for both the 'analysis' and the 'guess'.
%          The input datafile is the result of running 'obs_diag'.
%
% ddir     is an optional argument specifying the directory containing
%               the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
%
% ddir = 'plot';
% fit_ens_mean_time(ddir)

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

   % Be aware that the low order models override the varnames
   % and Regions variables defined above.

else
   error(sprintf('%s cannot be found.', datafile))
end

% Set up a structure with all the plotting components

plotdat.level    = level;
plotdat.toff     = toff;
plotdat.ylabel   = 'RMSE';
plotdat.nregions = length(Regions);
plotdat.nvars    = length(varnames);
plotdat.flavor   = 'Ens Mean';

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
         main = sprintf('%s %s',plotdat.flavor,string1);
      otherwise
         ges  = sprintf('%sges_times_%04dmb.dat',varnames{ivar},level);
         anl  = sprintf('%sanl_times_%04dmb.dat',varnames{ivar},level);
         main = sprintf('%s %s %d hPa',plotdat.flavor,string1,plotdat.level);
   end

   plotdat.ges     = ges;
   plotdat.anl     = anl;

   % plot each region

   figure(ivar); clf; 

   for iregion = 1:length(Regions),
      plotdat.title  = Regions{iregion};
      plotdat.region = iregion;
      myplot(plotdat);
   end

   CenterAnnotation(main);  % One title in the middle
   BottomAnnotation(ges);   % annotate filename at bottom

   % create a postscript file

   psfname = sprintf('%s_ens_mean_time.ps',plotdat.varname);
   print(ivar,'-dpsc',psfname);

end

path(startpath); % restore MATLABPATH to original setting

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------

function myplot(plotdat)
p1 = load(plotdat.ges); p = SqueezeMissing(p1); 
a1 = load(plotdat.anl); a = SqueezeMissing(a1); 

xp = p(:,1) + p(:,2)/86400 + plotdat.toff;
xa = a(:,1) + a(:,2)/86400 + plotdat.toff;

offset = 3;  % columns 1,2 are time, 3=mean, 4=spread, 5=numobs

count  = offset+(plotdat.region-1)*3;
yp     = p(:,count);
ya     = a(:,count);

gmean = mean(yp); gstring = sprintf('guess;    mean=%.3f',gmean);
amean = mean(ya); astring = sprintf('anaylsis; mean=%.3f',amean);

subplot(2,2,plotdat.region)
   plot(xp,yp,'k+-',xa,ya,'ro-','LineWidth',1.5)
   grid
   ax = axis; ax(3) = 0.0; axis(ax);
   datetick('x',1)
   ylabel(plotdat.ylabel, 'fontsize', 10);
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
