function fit_ens_mean_vertical(ddir)
% fit_ens_mean_vertical(ddir)
%
% ddir     is an optional argument specifying the directory containing
%               the data files as preprocessed by the support routines.
%
% USAGE: if the preprocessed data files are in a directory called 'plot'
%
%Wanl_ver_avedat ddir = 'plot';
% fit_ens_mean_vertical(ddir)

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

% This ensures the directory with the datafiles 
% is in Matlab's search path.

% Ensures the datafiles exist.
if ( nargin > 0 )
   datafile = fullfile(ddir,'ObsDiagAtts');
   TGuessFname = fullfile(ddir,'Tges_ver_ave.dat');
   TAnalyFname = fullfile(ddir,'Tanl_ver_ave.dat');
   WGuessFname = fullfile(ddir,'Wges_ver_ave.dat');
   WAnalyFname = fullfile(ddir,'Wanl_ver_ave.dat');
else
   datafile = 'ObsDiagAtts';
   TGuessFname = 'Tges_ver_ave.dat';
   TAnalyFname = 'Tanl_ver_ave.dat';
   WGuessFname = 'Wges_ver_ave.dat';
   WAnalyFname = 'Wanl_ver_ave.dat';
end
if ( exist(TGuessFname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.', TGuessFname))
end
if ( exist(WGuessFname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.', WGuessFname))
end
if ( exist(TAnalyFname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.', TAnalyFname))
end
if ( exist(WAnalyFname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.', WAnalyFname))
end

% subplot('position', [0.1,0.6,0.35,0.35])  ==? subplot(2,2,1)
% subplot('position', [0.6,0.6,0.35,0.35])  ==? subplot(2,2,2)
% subplot('position', [0.1,0.1,0.35,0.35])  ==? subplot(2,2,3)
% subplot('position', [0.6,0.1,0.35,0.35])  ==? subplot(2,2,4)

%----------------------------------------------------------------------
% Get attributes from obs_diag run.
%----------------------------------------------------------------------
  
if ( exist(datafile) == 2 )

   eval(datafile)

   temp = datenum(obs_year,obs_month,obs_day);
   toff = temp - round(t1); % determine temporal offset (calendar base)
   day1 = datestr(t1+toff,'yyyy-mm-dd HH');
   dayN = datestr(tN+toff,'yyyy-mm-dd HH');

else
   error(sprintf('%s cannot be found.', datafile))
end

main = sprintf('Ensemble Mean %s - %s',day1,dayN);

%----------------------------------------------------------------------
figure(1); clf; % Temperature
%----------------------------------------------------------------------

pv = load(TGuessFname); p_v = SqueezeMissing(pv); yp_v = p_v(:,1);
av = load(TAnalyFname); a_v = SqueezeMissing(av); ya_v = a_v(:,1);

ylab   = 'Pressure (hPa)';
xlab   = 'Temperature RMSE';

% Try to figure out intelligent axis limits
xdatarr = [p_v(:,2:2:8)  a_v(:,2:2:8)];      % concatenate all data
xlims   = [0.0 max(xdatarr(:))]; % limits of all data
ydatarr = [p_v(:,1) a_v(:,1)];               % concatenate all data
ylims   = [min(ydatarr(:)) max(ydatarr(:))]; % limits of all data
axlims  = [floor(xlims(1)) ceil(xlims(2)) round(ylims)];

region = 'Northern Hemisphere';
myplot(1, p_v(:,2), yp_v, a_v(:,2), ya_v, xlab, ylab, region, axlims)
region = 'Southern Hemisphere';
myplot(2, p_v(:,4), yp_v, a_v(:,4), ya_v, xlab, ylab, region, axlims)
region = 'Tropics';
myplot(3, p_v(:,6), yp_v, a_v(:,6), ya_v, xlab, ylab, region, axlims)
region = 'North America';
myplot(4, p_v(:,8), yp_v, a_v(:,8), ya_v, xlab, ylab, region, axlims)

CenterAnnotation(main)

%----------------------------------------------------------------------
figure(2); clf; % Windspeed
%----------------------------------------------------------------------

pv = load(WGuessFname); p_v = SqueezeMissing(pv); yp_v = p_v(:,1);
av = load(WAnalyFname); a_v = SqueezeMissing(av); ya_v = a_v(:,1);

ylab   = 'Pressure (hPa)';
xlab   = 'Windspeed RMSE';

% Try to figure out intelligent axis limits
xdatarr = [p_v(:,2:2:8)  a_v(:,2:2:8)];      % concatenate all data
xlims   = [0.0 max(xdatarr(:))]; % limits of all data
ydatarr = [p_v(:,1) a_v(:,1)];               % concatenate all data
ylims   = [min(ydatarr(:)) max(ydatarr(:))]; % limits of all data
axlims  = [floor(xlims(1)) ceil(xlims(2)) round(ylims)];

region = 'Northern Hemisphere';
myplot(1, p_v(:,2), yp_v, a_v(:,2), ya_v, xlab, ylab, region, axlims)
region = 'Southern Hemisphere';
myplot(2, p_v(:,4), yp_v, a_v(:,4), ya_v, xlab, ylab, region, axlims) 
region = 'Tropics';
myplot(3, p_v(:,6), yp_v, a_v(:,6), ya_v, xlab, ylab, region, axlims) 
region = 'North America';
myplot(4, p_v(:,8), yp_v, a_v(:,8), ya_v, xlab, ylab, region, axlims)

CenterAnnotation(main)

print -f1 -dpsc t_vertical.ps 
print -f2 -dpsc w_vertical.ps



function myplot(figpos,gx,gy,ax,ay,xlab,ylab,region,axlims)

subplot(2,2,figpos)
plot(gx,gy,'k+-',ax,ay,'ro-','LineWidth',1.5)
axis(axlims)
grid
set(gca,'YDir', 'reverse')
title({region}, 'FontSize', 14, 'FontWeight', 'bold' )
ylabel(ylab, 'fontsize', 10)
xlabel(xlab, 'fontsize', 10)
if   isempty(strfind(lower(xlab),'wind')) 
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



function CenterAnnotation(top)
subplot('position',[0.48 0.48 0.04 0.04])
axis off
h = text(0.5,0.5,top);
set(h,'HorizontalAlignment','center','VerticalAlignment','bottom',...
   'FontSize',12,'FontWeight','bold')

