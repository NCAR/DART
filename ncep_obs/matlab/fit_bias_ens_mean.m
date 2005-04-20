function fit_bias_ens_mean(ddir)
% fit_bias_ens_mean(ddir)
%
% ddir     is an optional argument specifying the directory containing
%               the data files as preprocessed by the support routines.
%
% USAGE:
%
% fit_bias_ens_mean('plot')
%
% Remember you can click and drag the legends ...

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

% Ensures the datafiles exist.
if ( nargin > 0 )
   TGuessFname = fullfile(ddir,'Tges_ver_ave_bias.dat');
   TAnalyFname = fullfile(ddir,'Tanl_ver_ave_bias.dat');
   WGuessFname = fullfile(ddir,'Wges_ver_ave_bias.dat');
   WAnalyFname = fullfile(ddir,'Wanl_ver_ave_bias.dat');
else
   TGuessFname = 'Tges_ver_ave_bias.dat';
   TAnalyFname = 'Tanl_ver_ave_bias.dat';
   WGuessFname = 'Wges_ver_ave_bias.dat';
   WAnalyFname = 'Wanl_ver_ave_bias.dat';
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

%----------------------------------------------------------------------
figure(1); clf; % Temperatures first, then
%----------------------------------------------------------------------

p_v  = load(TGuessFname);
a_v  = load(TAnalyFname);
yp_v = p_v(:,1);
ya_v = a_v(:,1);

ylab   = 'Pressure (hPa)';
xlab   = 'Bias';
top    = 'T Bias, Ensemble Mean';

% Try to figure out intelligent axis limits
xdatarr = [p_v(:,2:2:8)  a_v(:,2:2:8)];      % concatenate all data
xlims   = [min(xdatarr(:)) max(xdatarr(:))]; % limits of all data
ydatarr = [p_v(:,1) a_v(:,1)];               % concatenate all data
ylims   = [min(ydatarr(:)) max(ydatarr(:))]; % limits of all data
axlims  = [floor(xlims(1)) ceil(xlims(2)) round(ylims)];

region = 'Northern Hemisphere';
myplot(1, p_v(:,2), yp_v, a_v(:,2), ya_v, xlab, ylab, region, top, axlims)
region = 'Southern Hemisphere';
myplot(2, p_v(:,4), yp_v, a_v(:,4), ya_v, xlab, ylab, region, top, axlims)
region = 'Tropics';
myplot(3, p_v(:,6), yp_v, a_v(:,6), ya_v, xlab, ylab, region, top, axlims)
region = 'North America';
myplot(4, p_v(:,8), yp_v, a_v(:,8), ya_v, xlab, ylab, region, top, axlims)

%----------------------------------------------------------------------
figure(2); clf; % Winds
%----------------------------------------------------------------------

p_v  = load(WGuessFname);
a_v  = load(WAnalyFname);
yp_v = p_v(:,1);
ya_v = a_v(:,1);

ylab   = 'Pressure (hPa)';
xlab   = 'm/s';
top    = 'Windspeed Bias, Ensemble Mean';

% Try to figure out intelligent axis limits
xdatarr = [p_v(:,2:2:8)  a_v(:,2:2:8)];      % concatenate all data
xlims   = [min(xdatarr(:)) max(xdatarr(:))]; % limits of all data
ydatarr = [p_v(:,1) a_v(:,1)];               % concatenate all data
ylims   = [min(ydatarr(:)) max(ydatarr(:))]; % limits of all data
axlims  = [floor(xlims(1)) ceil(xlims(2)) round(ylims)];

region = 'Northern Hemisphere';
myplot(1, p_v(:,2), yp_v, a_v(:,2), ya_v, xlab, ylab, region, top, axlims)
region = 'Southern Hemisphere';
myplot(2, p_v(:,4), yp_v, a_v(:,4), ya_v, xlab, ylab, region, top, axlims)
region = 'Tropics';
myplot(3, p_v(:,6), yp_v, a_v(:,6), ya_v, xlab, ylab, region, top, axlims)
region = 'North America';
myplot(4, p_v(:,8), yp_v, a_v(:,8), ya_v, xlab, ylab, region, top, axlims)

%----------------------------------------------------------------------
print -f1 -dpsc t_bias.ps 
print -f2 -dpsc w_bias.ps
%----------------------------------------------------------------------

function myplot(figpos,gx,gy,ax,ay,xlab,ylab,region,top,axlims)

subplot(2,2,figpos)
   plot(gx,gy,'k+-',ax,ay,'ro-','LineWidth',2.0)
   axis(axlims)
   grid
   set(gca,'YDir', 'reverse')
   hold on; plot([0 0],[axlims(3) axlims(4)],'k-')
   title({top,region}, 'FontSize', 14, 'FontWeight', 'bold' )
   ylabel(ylab, 'fontsize', 10)
   xlabel(xlab, 'fontsize', 10)
   h = legend('guess', 'analysis','Location','best');
   legend(h,'boxoff')
