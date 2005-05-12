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
   datafile = fullfile(ddir,'Tanl_times_level');
   Tfname = fullfile(ddir,'Tges_ver_ave.dat');
   Wfname = fullfile(ddir,'Wges_ver_ave.dat');
else
   datafile = 'Tanl_times_level';
   Tfname = 'Tges_ver_ave.dat';
   Wfname = 'Wges_ver_ave.dat'; 
end

%----------------------------------------------------------------------
% Get attributes from obs_diag run.
%----------------------------------------------------------------------

if ( exist(datafile) == 2 )

   eval(datafile)

   temp   = datenum(obs_year,obs_month,obs_day);
   toff = temp - round(t1); % determine temporal offset (calendar base)
   day1 = datestr(t1+toff,'yyyy-mm-dd HH');
   dayN = datestr(tN+toff,'yyyy-mm-dd HH');

else
   error(sprintf('%s cannot be found.', datafile))
end

if ( exist(Tfname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',Tfname))
end
if ( exist(Wfname,'file') ~= 2 )
   error(sprintf('%s does not seem to exist.',Wfname))
end

linewidth = 2.0;

%----------------------------------------------------------------------
figure(1); clf;   % Temperature
%----------------------------------------------------------------------
switch obs_select
   case 1,
      string1 = sprintf('%s (all data)', 'T');
   case 2, 
      string1 = sprintf('%s (RaObs)', 'T');
   otherwise,
      string1 = sprintf('%s (ACARS,SATWND)', 'T');
end

pv      = load(Tfname); p_v = SqueezeMissing(pv);

yp_v     = p_v(:,1); 
nT_NH    = p_v(:,3);
nT_SH    = p_v(:,5);
nT_TR    = p_v(:,7);
nT_NA    = p_v(:,9);

subplot('position', [0.1,0.6,0.35,0.35])
plot(nT_NH, yp_v, 'gs-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Northern Hemisphere','fontsize', 12,'FontWeight','bold')
ylabel('Pressure (hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

subplot('position', [0.6,0.6,0.35,0.35])
plot(nT_SH, yp_v, 'bd-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Southern Hemisphere','fontsize', 12,'FontWeight','bold')
ylabel('Pressure (hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

subplot('position', [0.1,0.1,0.35,0.35])
plot(nT_TR, yp_v, 'ro-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Tropics','fontsize', 12,'FontWeight','bold')
ylabel('Pressure (hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

subplot('position', [0.6,0.1,0.35,0.35])
plot(nT_NA, yp_v, 'k+-', 'LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('North America','fontsize', 12,'FontWeight','bold')
ylabel('Pressure (hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

main = sprintf('# of %s %s -> %s',string1,day1,dayN);

CenterAnnotation(main)

% disp('Pausing, hit any key ...'); pause
print -dpsc t_num_vert.ps

%----------------------------------------------------------------------
figure(2); clf; % All regions on one large figure
%----------------------------------------------------------------------

h = plot(nT_NH, yp_v, 'gs-', ...
         nT_SH, yp_v, 'bd-', ...
         nT_TR, yp_v, 'ro-', ...
         nT_NA, yp_v, 'k+-', 'LineWidth', linewidth);
grid
set(gca,'YDir', 'reverse')
title(main,'FontSize',14,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

h = legend('Northern Hemisphere','Southern Hemisphere', ...
       'Tropics','North America','Location','SouthEast');
legend(h,'boxoff')

print -dpsc -append t_num_vert.ps

%----------------------------------------------------------------------
figure(3); clf; % Wind Observations ... Individual regions 
%----------------------------------------------------------------------
switch obs_select
   case 1,
      string1 = sprintf('%s (all data)', 'T');
   case 2, 
      string1 = sprintf('%s (RaObs)', 'T');
   otherwise,
      string1 = sprintf('%s (ACARS,SATWND)', 'T');
end

pv   = load(Wfname); p_v = SqueezeMissing(pv);

yp_v  = p_v(:,1);
nW_NH = p_v(:,3);
nW_SH = p_v(:,5);
nW_TR = p_v(:,7);
nW_NA = p_v(:,9);

subplot('position', [0.1,0.6,0.35,0.35])
plot(nW_NH ,yp_v,'gs-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Northern Hemisphere','fontsize', 12,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

subplot('position', [0.6,0.6,0.35,0.35])
plot(nW_SH ,yp_v,'bd-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Southern Hemisphere','fontsize', 12,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

subplot('position', [0.1,0.1,0.35,0.35])
plot(nW_TR ,yp_v,'ro-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('Tropics','fontsize', 12,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

subplot('position', [0.6,0.1,0.35,0.35])
plot(nW_NA ,yp_v,'k+-','LineWidth', linewidth)
grid
set(gca,'YDir', 'reverse')
title('North America','fontsize', 12,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

main = sprintf('# of %s %s -> %s',string1,day1,dayN);

CenterAnnotation(main)

print -dpsc w_num_vert.ps

%----------------------------------------------------------------------
figure(4); clf; % All regions on one large figure
%----------------------------------------------------------------------

h = plot(nT_NH, yp_v, 'gs-', ...
         nT_SH, yp_v, 'bd-', ...
         nT_TR, yp_v, 'ro-', ...
         nT_NA, yp_v, 'k+-', 'LineWidth', linewidth);
grid
set(gca,'YDir', 'reverse')
title(main,'FontSize',14,'FontWeight','bold')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('observation count','fontsize',10)

h = legend('Northern Hemisphere','Southern Hemisphere', ...
       'Tropics','North America','Location','SouthEast');
legend(h,'boxoff')

print -dpsc -append w_num_vert.ps


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
set(h,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
   'FontSize',12,'FontWeight','bold')

