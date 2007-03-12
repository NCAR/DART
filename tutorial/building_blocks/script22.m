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

axis([-4 4 0 0.8]);
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
text(-3.6, 0.75, '"TRUE" Prior PDF', 'fontsize', 24);

% Overlay the S.D. of distribution
sdx = [0, 0.6] - 2.0;
sdy = [0.2, 0.2];
%h_sd = plot(sdx, sdy, 'linewidth', 3);
%set(h_sd, 'color', [0 0.73 0]);
%h_sd_label = text(-2.6, 0.2, 'S.D.', 'fontsize', 24);

pause;
print -depsc s22f01.eps;

% Put on an ensemble shifted to the right by 2
x_prior_bad_mean = x_prior + 3.0;
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior_bad_mean = (x_prior_bad_mean - mean(x_prior_bad_mean)) * 0.6 / std(x_prior_bad_mean) + 1.0
h_plot = plot(x_prior_bad_mean, y_prior, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 1)
set(h_plot, 'color', [0 0.73 0]);
prior_bad_mean = normpdf(x, -2.0 + 3.0, 0.6);
h_prior = plot(x, prior_bad_mean, 'g--', 'linewidth', 3)
set(h_prior, 'color', [0 0.73 0]);
text(0.0, 0.75, 'Error in Mean (from model)', 'fontsize', 24);
% Overlay the S.D. of distribution
sdx = [0, 0.6] - 2.0 + 3.0;
sdy = [0.2, 0.2];
%h_sd = plot(sdx, sdy, 'g--', 'linewidth', 3);
%set(h_sd, 'color', [0 0.73 0]);
%h_sd_label = text(-2.6, 0.2, 'S.D.', 'fontsize', 24);

pause;
print -depsc s22f02.eps;

% Improve consistency by inflating
x_inf = sqrt(9) * (x_prior_bad_mean - mean(x_prior_bad_mean)) + mean(x_prior_bad_mean)
y_prior = y_prior + 0.05;
h_plot = plot(x_inf, y_prior, 'm*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
prior_inf = normpdf(x, -2.0 + 3.0, sqrt(9) * 0.6);
h_prior = plot(x, prior_inf, 'm', 'linewidth', 3)
text(1.5, 0.26, 'Variance Inflated', 'fontsize', 24);
% Overlay the S.D. of distribution
%sdx = [0, sqrt(9) * 0.6] - 2.0 + 3.0;
%sdy = [0.3, 0.3];
%h_sd = plot(sdx, sdy, 'm', 'linewidth', 3);

pause;
print -depsc s22f03.eps;
