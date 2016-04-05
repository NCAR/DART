% A time series from an L63 assim
% Identity Obs. every 24 hours with 4.0 SD
% 4x20 ensemble, prior plotted green, obs red
% Long background climatology in light blue 

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

% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 5.8 5.8]); 

load l63_attract;
grid on;
hold on;
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
%set(gcf, 'yticklabel', []);
%set(gcf, 'Xticklabel', []);
axis([-20, 20, -20, 20, 10, 40]);

% SHOULD USE 10 as upper bound, but crashes at CGD
for i = 1:10
   lower = (i - 1) * 300 + 1;
   upper = lower + 299;
   h1 = plot3(l63_attract(lower:upper, 1), l63_attract(lower:upper, 2), l63_attract(lower:upper, 3), '*', 'markersize', 2);
   i
   pause;
   %set(h1, 'markersize', 2)
end

% Load the ensemble
load l63_ens;
ens = l63_ens;

% Load the obs
load l63_obs;
obs = l63_obs;

for i = 2:10
   low = (i - 1) * 80 + 1;
   hi = low + 79;
   h = plot3(ens(low:hi, 2), ens(low:hi, 3), ens(low:hi, 4), 'g.');
   set(h, 'color', [0 0.73 0]);
   view(45, 45);
   set(h, 'markersize', 36);
   h_obs = plot3(obs(i:i, 1), obs(i:i, 2), obs(i:i, 3), 'r.');
   set(h_obs, 'markersize', 36);

   pause;
   % Setup the printing characteristics
   set(gcf, 'PaperPositionMode', 'auto');
   %set(gcf, 'PaperOrientation', 'landscape');
   set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
   outname = ['s19f0', num2str(i)]
   print(gcf, '-depsc', outname);

   % Now show the inflated ensemble and the ensemble mean in some other color
   lamda = 1.8
   for j = 2:4
      ens_mean(j) = mean(ens(low:hi, j));
      inf_ens(1:80, j) = lamda * (ens(low:hi, j) - ens_mean(j)) + ens_mean(j);
   end
   h_inf = plot3(inf_ens(:, 2), inf_ens(:, 3), inf_ens(:, 4), 'm.');
   set(h_inf, 'markersize', 36);
   h_mean = plot3(ens_mean(2), ens_mean(3), ens_mean(4), 'm*');
   set(h_mean, 'markersize', 36);

   pause;
   % Would need to print this set, too
   outname = ['s19af0', num2str(i)]
   print(gcf, '-depsc', outname);

   set(h, 'visible', 'off');
   set(h_obs, 'visible', 'off');
   set(h_inf, 'visible', 'off');
   set(h_mean, 'visible', 'off');
end
