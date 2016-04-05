%

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

% Lots of interesting non-Gaussian stuff can happen too
x = -5:0.01:5;
prior = (normpdf(x, -1.5) + normpdf(x, 1.5)) ./ 2.0;
obs = normpdf(x, 1.0, 2.0);
product = prior .* obs;
hhh = plot(x, prior, 'g', 'linewidth', 3)
set(hhh, 'color', [0, 0.73, 0]);


axis([-4 4 0 0.3]);
set(gca, 'fontsize', 24);
text(-3.4, 0.2, 'Prior PDF', 'fontsize', 24);
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
% PRINTING IS SCREWED UP WITH CURRENT MATLAB VERSION, DO MANUALLY
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print -depsc s03f01.eps;

plot(x, obs, 'r', 'linewidth', 3)
text(-1.6, 0.075, 'Obs. Likelihood', 'fontsize', 24);
pause;
print -depsc s03f02.eps;

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
text(-1.1, 0.27, 'Posterior PDF', 'fontsize', 24);
print -depsc s03f03.eps;

