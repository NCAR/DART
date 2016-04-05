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
set(gcf, 'position', [1 1 6.5 6.5]);

x = d(1:21, 1);
y10 = d(1:21, 2);
y20 = d(22:42, 2);
y40 = d(43:63, 2);
y80 = d(64:84, 2);

set(gca, 'Position', [0.18 0.14 0.7550 0.7550]);
plot(x, y10, 'b', 'linewidth', 3.0);
hold on;
set(gca, 'linewidth', 3);
set(gca, 'fontsize', 24);
xlabel('True Correlation');
ylabel('Expected |Sample Correlation|');
hhh = plot(x, y20, 'g', 'linewidth', 3.0);
set(hhh, 'color', [0 0.73 0]);
plot(x, y40, 'r', 'linewidth', 3.0);
plot(x, y80, 'm', 'linewidth', 3.0);
grid on;
legend('10 Members', '20 Members', '40 Members', '80 Members', 4);
plot(x, x, 'k', 'linewidth', 3.0);
% Replot these to make them overlay in desired order after setting legend
plot(x, y80, 'm', 'linewidth', 3.0);
plot(x, y40, 'r', 'linewidth', 3.0);
hhh = plot(x, y20, 'g', 'linewidth', 3.0);
set(hhh, 'color', [0 0.73 0]);
plot(x, y10, 'b', 'linewidth', 3.0);

pause
% Setup the printing characteristics                                                     
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s15f01.eps

axis([0 0.5 0 0.5]);

print -depsc s15f02.eps

