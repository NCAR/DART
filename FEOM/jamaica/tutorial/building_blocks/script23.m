% Looking at inflation to deal with errors

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

% Begin by setting random seed to a nice controlled
% initial value that gives nice plots
randn('state', 0);

x_prior = [-2.2 -1.3 -0.4 1.2 2.2];
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior = (x_prior - mean(x_prior)) * 0.6 / std(x_prior) - 2.0

% Plot the prior distribution
y_prior = [0.02, 0.02, 0.02, 0.02, 0.02];
h_plot = plot(x_prior, y_prior, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
set(h_plot, 'color', [0 0.73 0]);

axis([-4 0 0 1.4]);
set(gca, 'fontsize', 24);
p_text = text(-2.0, 0.1, 'Prior Ensemble', 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');
hold on;

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');

% Plot a gaussian fit to ensemble
set(p_text, 'visible', 'off');
x = -5:0.01:5;
prior = normpdf(x, -2.0, 0.6);
h_prior = plot(x, prior, 'g', 'linewidth', 3)
set(h_prior, 'color', [0 0.73 0]);
text(-3.7, 0.75, '"TRUE" Prior PDF', 'fontsize', 24);

pause;
print -depsc s23f01.eps;

% Put on an ensemble with insufficient variance
x_prior_def = 0.5 * (x_prior - mean(x_prior)) + mean(x_prior);
y_prior_def = y_prior + 0.25;
h_plot = plot(x_prior_def, y_prior_def, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 1)
set(h_plot, 'color', [0 0.73 0]);

% Plot a gaussian fit to ensemble
set(p_text, 'visible', 'off');
x = -5:0.01:5;
prior = normpdf(x, -2.0,  0.5 * 0.6);
h_prior = plot(x, prior, 'g--', 'linewidth', 3)
set(h_prior, 'color', [0 0.73 0]);
text(-3.8, 1.2, 'Variance Deficient PDF', 'fontsize', 24);

pause;
print -depsc s23f02.eps;


