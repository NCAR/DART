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
%h1 = plot3(l63_attract(1:3000, 1), l63_attract(1:3000, 2), l63_attract(1:3000, 3), '*');
h1 = plot3(l63_attract(1:1000, 1), l63_attract(1:1000, 2), l63_attract(1:1000, 3), '*');
set(h1, 'markersize', 2)
grid on;
hold on;
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
%set(gcf, 'yticklabel', []);
%set(gcf, 'Xticklabel', []);
axis([-15, 15, -20, 20, 10, 40]);

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
   outname = ['s17f0', num2str(i)]
   print(gcf, '-depsc', outname);

   set(h, 'visible', 'off');
   set(h_obs, 'visible', 'off');
end
