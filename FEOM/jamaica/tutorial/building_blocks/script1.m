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

% First 1-D case, just two funky distributions slightly skewed, take product
x = -6.0:0.01:6.0;
% Make a somewhat skewed prior distribution
prior = (normpdf(x, -2, 1.0) + normpdf(x, 0.5, 1.5)) ./ 2;
hhh = plot(x, prior, 'g', 'linewidth', 3)
set(hhh, 'color', [0 0.73 0]);
axis([-6 6 0 0.25]);
set(gca, 'fontsize', 24);
text(-4.7, 0.2, 'Prior PDF', 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');
pause;

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s01f01.eps;

hold on;
obs = (normpdf(x, 2, 1.5) + normpdf(x, 0, 2.0)) ./ 2;
plot(x, obs, 'r', 'linewidth', 3)
text(2.6, 0.18, 'Obs. Likelihood', 'fontsize', 24);
pause;
print -depsc s01f02.eps;

product = prior .* obs;
plot(x, product, 'b--', 'linewidth', 3);
ph = text(-2.1, 0.04, 'Product (Numerator)', 'fontsize', 24);
pause;
print -depsc s01f03.eps;

h = area(x, product)
set(h, 'edgecolor', 'b')
set(h, 'facecolor', [0.0 0.75 1.0])
set(ph, 'visible', 'off');
text(-2.1, 0.04, 'Normalization (Denom.)', 'fontsize', 24);
pause;
print -depsc s01f04.eps;

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
text(1.30, 0.225, 'Posterior', 'fontsize', 24);

print -depsc s01f05.eps;
