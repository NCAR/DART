% Plot errors associated with finite size correlation
% correl_error.f90 in system_simulation in DART generates input

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

load correl_errors_stripped;
d = correl_errors_stripped;

% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 9.0 4.5]);

% Compute a Gaspari Cohn localization
x = -2000:100:2000;
half_width = 1000;

z = abs(x);
y(1:size(z, 2)) = 0.0;
for i = 1:size(z, 2)
   if (z(i) <= 2.0 * half_width && z(i) >= half_width) 
      r = z(i) / half_width;
      y(i) = r^5 / 12 - r^4 / 2 + r^3 * 5/8 + r^2 * 5/3 - 5*r + 4 - (half_width * 2) / (3*z(i));
   end
   if (z(i) < half_width)
      r = z(i) / half_width;
      y(i) = r^5 * (-1/4) + r^4 / 2 + r^3 * 5/8 - r^2 * 5/3 + 1;
   end
end

h_gc = plot(x, y, 'linewidth', 3.0)
axis([-2000, 2000, 0, 1]);
hold on;
set(gca, 'Position', [0.18 0.24 0.7550 0.6850]);
set(gca, 'linewidth', 3);
set(gca, 'fontsize', 24);
xlabel('Distance from Observation ');
ylabel('Regression Weight');
grid on;

pause

% Setup the printing characteristics                                                     
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s16f01.eps

figure(2)
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 9.0 3.0]);
hold on;
h_gc = plot(x, y, 'linewidth', 3.0);
grid on;
axis([-2100, 2100, 0, 1.05]);
set(gca, 'Position', [0.18 0.24 0.7550 0.6850]);
set(gca, 'linewidth', 3);
set(gca, 'fontsize', 24);
xlabel('Distance from Observation ');
ylabel('Weight');

% Setup the printing characteristics                                                     
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s16f02.eps

pause

set(h_gc, 'visible', 'off')
xb = [-2000 -2000 2000 2000];
yb = [0 .99 .99 0];
h_b = plot(xb, yb, 'linewidth', 3.0)

print -depsc s16f03.eps

pause

set(h_b, 'visible', 'off')
xrb = [-2000 -1000 1000 2000];
yrb = [0 1 1 0];
plot(xrb, yrb, 'linewidth', 3.0)

print -depsc s16f04.eps

