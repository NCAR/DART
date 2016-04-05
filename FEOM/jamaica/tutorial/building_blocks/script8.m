% Ensemble adjustment filter to maintain some structure (this is the one I use most)

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

x_prior = [-2.2 -1.8 -1.4 1.8 2.2];
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior = (x_prior - mean(x_prior)) * 1.2 / std(x_prior) - 1.0

% Plot the prior distribution
y_prior = [0.05, 0.05, 0.05, 0.05, 0.05];
h_plot = plot(x_prior, y_prior, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
set(h_plot, 'color', [0 0.73 0]);

axis([-4 4 0 0.6]);
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
pause

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s08f01.eps;

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
product = prior .* obs;
h_prior_plot = plot(x, prior, 'g', 'linewidth', 3)
set(h_prior_plot, 'color', [0 0.73 0]);
h_prior_text = text(-3.0, 0.3, 'Prior PDF', 'fontsize', 24);
pause;
print -depsc s08f02.eps;

h_obs_plot = plot(x, obs, 'r', 'linewidth', 3)
h_obs_text = text(1.6, 0.4, 'Obs. Likelihood', 'fontsize', 24);
pause;
print -depsc s08f03.eps;

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
text(-1.9, 0.55, 'Posterior PDF', 'fontsize', 24);
pause
print -depsc s08f04.eps;

% Now remove the obs. likelihood and the prior curve to keep things clean
set(h_obs_plot, 'visible', 'off');
set(h_prior_plot, 'visible', 'off');
set(p_text, 'visible', 'off');
set(h_prior_text, 'visible', 'off');
set(h_obs_text, 'visible', 'off');
pause;
print -depsc s08f05.eps;

% First, need to compute the mean and variance of the blue
posterior_var = 1.0 / (1.0 / (1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_mean = posterior_var * (-1.0 / (1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_sd = sqrt(posterior_var);

% First shift the mean to have posterior mean
x_posterior = x_prior - mean(x_prior) + posterior_mean
y_posterior = y_prior + 0.1;
h_plot = plot(x_posterior, y_posterior, 'b*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
text(-2.5, 0.2, 'Mean Shifted', 'fontsize', 24);

% plot update lines between each pair
for i = 1:5
   xu = [x_prior(i) x_posterior(i)];
   yu = [y_prior(i) y_posterior(i)];
   plot(xu, yu)
end

pause;
print -depsc s08f06.eps;

% Now adjust the variance by linear contraction
x_final = (x_posterior - posterior_mean) * sqrt(posterior_var / var(x_posterior)) + posterior_mean
y_final = y_posterior + 0.1;
h_plot = plot(x_final, y_final, 'b*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
text(-2.0, 0.3, 'Variance Adjusted', 'fontsize', 24);

% plot update lines between each pair
for i = 1:5
   xu = [x_final(i) x_posterior(i)];
   yu = [y_final(i) y_posterior(i)];
   plot(xu, yu)
end

print -depsc s08f07.eps;
