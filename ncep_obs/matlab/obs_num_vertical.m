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

% This ensures the datafiles exist. 
if ( nargin > 0 )
   Tfname = fullfile(ddir,'Tges_ver_ave.dat');
   Wfname = fullfile(ddir,'Wges_ver_ave.dat');
else
   Tfname = 'Tges_ver_ave.dat';
   Wfname = 'Wges_ver_ave.dat'; 
end
if ( exist(Tfname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',Tfname))
end
if ( exist(Wfname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',Wfname))
end

% obsfit.m plot of the error of analysis and guess
figure(1); clf
p_v      = load(Tfname);

yp_v     = p_v(:,1); 
nT_NH    = p_v(:,3);
nT_SH    = p_v(:,5);
nT_TR    = p_v(:,7);
nT_NA    = p_v(:,9);

linewidth = 2.0;

subplot('position', [0.1,0.6,0.35,0.35])
plot(nT_NH, yp_v, 'gs-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Northern Hemisphere')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Temperature Observations','fontsize',10)

subplot('position', [0.6,0.6,0.35,0.35])
plot(nT_SH, yp_v, 'bd-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Southern Hemisphere')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Temperature Observations','fontsize',10)

subplot('position', [0.1,0.1,0.35,0.35])
plot(nT_TR, yp_v, 'ro-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Tropics')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Temperature Observations','fontsize',10)

subplot('position', [0.6,0.1,0.35,0.35])
plot(nT_NA, yp_v, 'k+-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('North America')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Temperature Observations','fontsize',10)

% disp('Pausing, hit any key ...'); pause
print -dpsc t_num_vert.ps

%----------------------------------------------------------------------
% All regions on one large figure
%----------------------------------------------------------------------

figure(2); clf

h = plot(nT_NH, yp_v, 'gs-', ...
         nT_SH, yp_v, 'bd-', ...
         nT_TR, yp_v, 'ro-', ...
         nT_NA, yp_v, 'k+-', 'LineWidth', linewidth);
grid
set(gca,'YDir', 'reverse')
title('# of Temperature Observations by Region','FontSize',14,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Observations','fontsize',10)

legend('Northern Hemisphere','Southern Hemisphere', ...
       'Tropics','North America','Location','SouthEast');

print -dpsc -append t_num_vert.ps

%----------------------------------------------------------------------
% Wind Observations ... Individual regions 
%----------------------------------------------------------------------

figure(3); clf
p_v   = load(Wfname);
yp_v  = p_v(:,1);
nW_NH = p_v(:,3);
nW_SH = p_v(:,5);
nW_TR = p_v(:,7);
nW_NA = p_v(:,9);

subplot('position', [0.1,0.6,0.35,0.35])
plot(nW_NH ,yp_v,'gs-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Northern Hemisphere')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('NH', 'fontsize', 10)
xlabel('# of Wind Observations','fontsize',10)


subplot('position', [0.6,0.6,0.35,0.35])
plot(nW_SH ,yp_v,'bd-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Southern Hemisphere')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Wind Observations','fontsize',10)

subplot('position', [0.1,0.1,0.35,0.35])
plot(nW_TR ,yp_v,'ro-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Tropics')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Wind Observations','fontsize',10)

subplot('position', [0.6,0.1,0.35,0.35])
plot(nW_NA ,yp_v,'k+-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('North America')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Wind Observations','fontsize',10)

print -dpsc w_num_vert.ps

%----------------------------------------------------------------------
% All regions on one large figure
%----------------------------------------------------------------------

figure(4); clf

h = plot(nT_NH, yp_v, 'gs-', ...
         nT_SH, yp_v, 'bd-', ...
         nT_TR, yp_v, 'ro-', ...
         nT_NA, yp_v, 'k+-', 'LineWidth', linewidth);
grid
set(gca,'YDir', 'reverse')
title('# of Wind Observations by Region','FontSize',14,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('# of Observations','fontsize',10)

legend('Northern Hemisphere','Southern Hemisphere', ...
       'Tropics','North America','Location','SouthEast');

print -dpsc -append w_num_vert.ps
