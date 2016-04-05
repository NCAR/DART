%% DART:script2

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Things are natural for Gaussians
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
product = prior .* obs;
hhh = plot(x, prior, 'g', 'linewidth', 3)
set(hhh, 'color', [0, 0.73, 0]);

axis([-4 4 0 0.6]);
set(gca, 'fontsize', 24);
text(-3.0, 0.3, 'Prior PDF', 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');

hold on;
pause;

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s02f01.eps;
 
plot(x, obs, 'r', 'linewidth', 3)
text(1.6, 0.4, 'Obs. Likelihood', 'fontsize', 24);
pause;
print -depsc s02f02.eps;

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
text(-1.9, 0.55, 'Posterior PDF', 'fontsize', 24);
print -depsc s02f03.eps;

