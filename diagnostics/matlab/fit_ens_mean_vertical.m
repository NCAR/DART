function fit_ens_mean_vertical(ddir)
% fit_ens_mean_vertical(ddir)
%
% Plots the RMS of the ensemble mean as a function of height for 
% several regions. The ensemble mean is averaged over a time period. 
% The calculations are done by 'obs_diag' - which generates data files 
% that are used by this plotting routine.
%
% the input data files are of the form *ges_ver_ave.dat,
% where the 'ave' refers to averaging over time. The first part of
% the file name is the name of the variable contained in the file.
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
% fit_ens_mean_vertical(ddir)

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
ptypes   = {'gs-','bd-','ro-','k+-'};    % for each region

%----------------------------------------------------------------------
% Get plotting metadata from obs_diag run.
%----------------------------------------------------------------------

if ( exist(datafile) == 2 )

   eval(datafile)

   if ( exist('plevel','var') == 0 )
      disp(sprintf('%s does not have multiple levels.', datafile))
      disp('It cannot be plotted with fit_ens_mean_vertical.')
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
plotdat.xlabel    = 'RMSE';
plotdat.linewidth = 2.0;

main = sprintf('Ensemble Mean %s - %s',plotdat.day1,plotdat.dayN);

%----------------------------------------------------------------------
% Loop around observation types
%----------------------------------------------------------------------

for ivar = 1:length(All_Level_Varnames),

   % set up a structure with all the plotting components

   plotdat.varname = All_Level_Varnames{ivar};

   switch obs_select
      case 1,
         string1 = sprintf('%s Ens Mean (all data)',     plotdat.varname);
      case 2,
         string1 = sprintf('%s Ens Mean (RaObs)',        plotdat.varname);
      otherwise,
         string1 = sprintf('%s Ens Mean (ACARS,SATWND)', plotdat.varname);
   end

   plotdat.ges  = sprintf('%s_ges_ver_ave.dat',All_Level_Varnames{ivar});
   plotdat.anl  = sprintf('%s_anl_ver_ave.dat',All_Level_Varnames{ivar});
   plotdat.main = sprintf('%s %sZ -- %sZ',string1,plotdat.day1,plotdat.dayN);
   plotdat      = SetLevels(plotdat);

   % plot by region

   figure(ivar); orient landscape; clf; wysiwyg

   for iregion = 1:length(Regions),
      plotdat.title  = Regions{iregion};
      plotdat.region = iregion;
      myplot(plotdat);
   end

   CenterAnnotation(plotdat.main)
   BottomAnnotation(plotdat.ges)

   % create a postscript file

   psfname = sprintf('%s_vertical.ps',plotdat.varname);
   print(ivar,'-dpsc',psfname);

end

path(startpath); % restore MATLABPATH to original setting

%----------------------------------------------------------------------
% 'Helper' functions
%----------------------------------------------------------------------

function myplot(plotdat)
regionindex = 2 + 2*(plotdat.region - 1);
pv     = load(plotdat.ges);
av     = load(plotdat.anl);
p_v    = SqueezeMissing(pv);
a_v    = SqueezeMissing(av);
guessX = p_v(:,regionindex);
analyX = a_v(:,regionindex);
levels = plotdat.levels;

% Try to figure out intelligent axis limits
indmax  = size(p_v,2);
xdatarr = [p_v(:,2:2:indmax)  a_v(:,2:2:indmax)];  % concatenate all data
xlims   = [0.0 max(xdatarr(:))];                   % limits of all data
axlims  = [floor(xlims(1)) ceil(xlims(2)) min(plotdat.ylims) max(plotdat.ylims)];

% sometimes there is no valid data, must patch axis limits
if (~isfinite(axlims(2)))
   axlims(2) =  1;
end

subplot(2,2,plotdat.region)
   plot(guessX,levels,'k+-',analyX,levels,'ro-','LineWidth',plotdat.linewidth)
   axis(axlims)
   grid
   set(gca,'YDir', plotdat.ydir)
   title( plotdat.title, 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold' )
   ylabel(plotdat.ylabel, 'FontSize', 10)
   xlabel(plotdat.xlabel, 'FontSize', 10)
   if   isempty(strfind(lower(plotdat.varname),'w')) 
      h = legend('guess', 'analysis','Location','East');
   else
      h = legend('guess', 'analysis','Location','SouthEast');
   end
   legend(h,'boxoff');



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

if ( exist(plotdat.ges,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',plotdat.ges))
end
if ( exist(plotdat.anl,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',plotdat.anl))
end

datmat = load(plotdat.ges);
guess  = datmat(:,1);  

datmat = load(plotdat.anl);
anl    = datmat(:,1);  

if (sum(abs(guess-anl)) > 0) 
   str1 = sprintf('levels in %s and %s do not match',plotdat.ges, plotdat.anl)
   error(str1)
end

levels = guess; clear guess anl 

plotstruct         = plotdat;
plotstruct.levels  = levels;
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

