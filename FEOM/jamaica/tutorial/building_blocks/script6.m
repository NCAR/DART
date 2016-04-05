% Generating a posterior sample 2: deterministic sampling

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

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
product = prior .* obs;

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
text(-1.8, 0.55, 'Posterior PDF', 'fontsize', 24);

axis([-3 4 0 0.6]);
set(gca, 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');
hold on;
pause

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s06f01.eps;

% Need to get a posterior ensemble
% Option 1, one example of exact mean and variance, kurtosis 3
% First, need to compute the mean and variance of the blue
posterior_var = 1.0 / (1.0 / (1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_mean = posterior_var * (-1.0 / (1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_sd = sqrt(posterior_var);

% Here's spacing for kurtosis 3, variance 1.0
x_posterior = [-1.956925 -0.7647584 -0.2830350 -7.387926E-02 -4.842975E-03 4.842975E-03 7.387926E-02 0.2830350 0.7647584 1.956925];
y_posterior = [0.05, 0.05, 0.05, 0.05 0.05 0.05 0.05 0.05 0.05 0.05];
% Adjust to mean and variance of posterior; kurtosis comes along for the ride
x_posterior = (x_posterior - mean(x_posterior)) * sqrt(var(x_posterior) / posterior_var) + posterior_mean
h_plot = plot(x_posterior, y_posterior, 'b*');
set(h_plot, 'markersize', 18);
set(h_plot, 'linewidth', 3);
h_text = text(-2.2, 0.10, 'Kurtosis 3', 'fontsize', 24);
pause
print -depsc s06f02.eps;


% Here's spacing for kurtosis 2, variance 1.0
x_posterior = [-1.716355 -1.040813 -0.6131358 -0.2999267 -7.030749E-02 7.030749E-02 0.2999267 0.6131358 1.040813 1.716355];
y_posterior = y_posterior + 0.12
% Adjust to mean and variance of posterior; kurtosis comes along for the ride
x_posterior = (x_posterior - mean(x_posterior)) * sqrt(var(x_posterior) / posterior_var) + posterior_mean
h_plot = plot(x_posterior, y_posterior, 'b*');
set(h_plot, 'markersize', 18);
set(h_plot, 'linewidth', 3);
h_text = text(-2.2, 0.22, 'Kurtosis 2', 'fontsize', 24);


print -depsc s06f03.eps;

