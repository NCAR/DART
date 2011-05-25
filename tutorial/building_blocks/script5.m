%% DART:script5  Generating a posterior sample 1: straight sampling of the continuous posterior distribution

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Begin by setting random seed to a nice controlled
% initial value that gives nice plots
randn('state', 0);

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
product = prior .* obs;

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
text(-1.2, 0.55, 'Posterior PDF', 'fontsize', 24);

axis([-2 3 0 0.6]);
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
% PRINTING IS SCREWED UP WITH CURRENT MATLAB VERSION, DO MANUALLY
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print -depsc s05f01.eps;

% Need to get a posterior ensemble
% Option 1, just do a random sample of the blue distribution
% First, need to compute the mean and variance of the blue
posterior_var = 1.0 / (1.0 / (1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_mean = posterior_var * (-1.0 / (1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_sd = sqrt(posterior_var);
x_posterior = normrnd(posterior_mean, posterior_sd, 5, 1);
y_posterior = [0.05, 0.05, 0.05, 0.05, 0.05];
h_plot = plot(x_posterior, y_posterior, 'b*');
set(h_plot, 'markersize', 18);
set(h_plot, 'linewidth', 3);
h_text = text(-1.5, 0.10, 'Random Sample', 'fontsize', 24);
pause
print -depsc s05f02.eps;

set(h_text, 'visible', 'off');
% Can adjust to have exact mean
x_mn = x_posterior - mean(x_posterior) + posterior_mean
y_mn = y_posterior + 0.08
h_plot = plot(x_mn, y_mn, 'b*');
set(h_plot, 'markersize', 18);
set(h_plot, 'linewidth', 3);
h_text = text(-1.8, 0.18, 'Random Sample; Exact Mean', 'fontsize', 24);
% plot update lines between each pair
for i = 1:5
   xu = [x_posterior(i) x_mn(i)];
   yu = [y_posterior(i) y_mn(i)];
   plot(xu, yu)
end
pause
print -depsc s05f03.eps;

set(h_text, 'visible', 'off');
% Can adjust to have exact mean and variance
x_var = (x_mn - posterior_mean) * sqrt(posterior_var / var(x_mn)) + posterior_mean
y_var = y_mn + 0.08
h_plot = plot(x_var, y_var, 'b*');
set(h_plot, 'markersize', 18);
set(h_plot, 'linewidth', 3);
h_text = text(-1.8, 0.26, 'Random Sample; Exact Mean and Var.', 'fontsize', 24);
% plot update lines between each pair
for i = 1:5
   xu = [x_mn(i) x_var(i)];
   yu = [y_mn(i) y_var(i)];
   plot(xu, yu)
end

print -depsc s05f04.eps;



