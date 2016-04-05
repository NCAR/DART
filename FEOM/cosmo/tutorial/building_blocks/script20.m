%% DART:script20 Illustrate the Bayesian update for the obs. space inflation param

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

% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [3 0.0 9.0 6.0]);
% Set the shape of the plot box
%pbaspect([2.4 1 1]);

% Take a look at using multiple axes here
r2 = subplot(2, 1, 2);
r1 = subplot(2, 1, 1);

% Set the prior y mean and sd and obs y and sd
y_prior_mean = 1.2;
y_prior_sd = 0.6;
y_obs = 2.0;
y_obs_sd = 0.6;


subplot(r1);
hold on;
grid on;
% Plot the observation distribution, the initial ensemble and fit
x_prior = [-0.2 0.7 1.6 3.2 4.2];
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior = (x_prior - mean(x_prior)) * y_prior_sd / std(x_prior) + y_prior_mean

% Plot the prior distribution
y_prior = [0.02, 0.02, 0.02, 0.02, 0.02];
h_plot = plot(x_prior, y_prior, 'g*')
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 20);
xlabel('Observation: y');
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
set(h_plot, 'color', [0 0.73 0]);

axis([-1 4 0 0.8]);
set(gca, 'fontsize', 20);

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, y_prior_mean, y_prior_sd);
obs = normpdf(x, y_obs, y_obs_sd);
h_prior = plot(x, prior, 'g--', 'linewidth', 3)
set(h_prior, 'color', [0 0.73 0]);
%text(-3.6, 0.7, 'Prior PDF', 'fontsize', 20);

% Overlay the S.D. of distribution
%sdx = [0, y_prior_sd] + y_prior_mean;
%sdy = [0.2, 0.2];
%h_sd = plot(sdx, sdy, 'linewidth', 3);
%set(h_sd, 'color', [0 0.73 0]);
%h_sd_label = text(-2.6, 0.2, 'S.D.', 'fontsize', 20);

plot(x, obs, 'r', 'linewidth', 3)
text(1.8, 0.72, 'Obs. Likelihood', 'fontsize', 20);
h_prior_text = text(0.1, 0.72, 'Prior PDF', 'fontsize', 20);


% Overlay the S.D. of distribution
%sdx = [-1 * y_obs_sd, 0] + y_obs;
%sdy = [0.2, 0.2];
%hhh = plot(sdx, sdy, 'r', 'linewidth', 3);
%text(2.1, 0.2, 'S.D.', 'fontsize', 20);

% Also have a prior distribution for lambda, normal here
subplot(r2);
hold on;
grid on;

% Position the marginal and joint plot boxes
%set(r1, 'Position', [0.13, 0.6111, 0.7750, 0.31]);
%set(r2, 'Position', [0.13, 0.25, 0.7750, 0.31]);

lam_prior_mean = 1.5;
lam_prior_sd = 0.25;
lam_prior = normpdf(x, lam_prior_mean, lam_prior_sd);
h_lam_prior = plot(x, lam_prior, 'g', 'linewidth', 3)
hold on;
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 20);
xlabel('Obs. Space Inflation Factor: \lambda');
set(h_lam_prior, 'color', [0 0.73 0]);
text(1.5, 1.8, 'Prior \lambda PDF', 'fontsize', 20);
axis([0 6 0 2.0]);

pause
% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s20f01.eps;


% Turn off some text
set(h_prior_text, 'visible', 'off');

% Compute actual distance between ensemble mean and observation
actual_dist = abs(y_prior_mean - y_obs);

% Pick an initial value of lambda to demonstrate the product
% Need to compute likelihood that yobs was observed given a value of lambda

for i = 1:3
   % Try starting with mean value from prior
   lambda = lam_prior_mean + (i - 2) * 3 * lam_prior_sd;

   % Show inflated values on top axis
   subplot(r1);
   x_prior_inf = sqrt(lambda) * (x_prior - mean(x_prior)) + mean(x_prior);
   y_prior_inf = y_prior + 0.1
   h_plot_inf = plot(x_prior_inf, y_prior_inf, 'g*');
   set(h_plot_inf, 'markersize', 18)
   set(h_plot_inf, 'linewidth', 3)
   set(h_plot_inf, 'color', [0 0.73 0]);
   % Also plot the expanded continuous distribution
   prior_inf = normpdf(x, y_prior_mean, y_prior_sd * sqrt(lambda));
   h_prior_inf(i) = plot(x, prior_inf, 'g', 'linewidth', 3)
   set(h_prior_inf(i), 'color', [0 0.73 0]);

   % Label the inflated pdf
   h_text_inf(i) = text(-0.95, 0.60, ['Inflated Prior \lambda = ', num2str(lambda)], 'fontsize', 18);

   % Compute the sd of the distance (expected value is 0 for unbiases)
   subplot(r2);
   dist_sd = sqrt(lambda * y_prior_sd^2 + y_obs_sd^2);
   % Probability of actual dist from normal is
   pr_this_dist = (1.0 / sqrt(dist_sd * pi)) * exp(-1.0 * actual_dist^2 / (2 * dist_sd^2));
   h_lam_obs = plot(lambda, pr_this_dist, 'r*', 'markersize', 20, 'linewidth', 3);
   % Compute prior value for this lambda for product
   prior_lam_val = normpdf(lambda, lam_prior_mean, lam_prior_sd);
   h_lam_post = plot(lambda, prior_lam_val * pr_this_dist, 'b*', 'markersize', 20, 'linewidth', 3);
   pause
   outname = ['s20f0' num2str(i+1)];
   print(gcf, '-depsc', outname);

   % Turn off the inflated ensemble and density
   set(h_prior_inf(i), 'visible', 'off');
   set(h_plot_inf, 'visible', 'off');
   set(h_text_inf(i) , 'visible', 'off');
   % Make the asterisk for the likelihood smaller
   set(h_lam_obs, 'markersize', 6);
   set(h_lam_post, 'markersize', 6);
end



pause;
print -depsc s20f05.eps;

% Loop to print additional range of values of likelihood and product
num_points = 1000;
prob_this_dist(1:1000) = 0.0;
normalization = 0.0;
for i = 1:num_points;
   % Compute the obs. likelihood for a set of values of lambda
   lambda(i) = i * 8 / num_points;
   % Compute the sd of the distance (expected value is 0 for unbiases)
   dist_sd = sqrt(lambda(i) * y_prior_sd^2 + y_obs_sd^2);
   % Probability of actual dist from normal is
   prob_this_dist(i) = (1.0 / sqrt(dist_sd * pi)) * exp(-1.0 * actual_dist^2 / (2 * dist_sd^2));
   % Compute prior value for this lambda for product
   prior_lam_val(i) = normpdf(lambda(i), lam_prior_mean, lam_prior_sd);
   post_lam_val(i) = prior_lam_val(i) * prob_this_dist(i);
   normalization = normalization + (8.0 / num_points) * post_lam_val(i);
end
plot(lambda, prob_this_dist, 'r', 'linewidth', 3);
plot(lambda, post_lam_val, 'b', 'linewidth', 3);
%plot(lambda, post_lam_val / normalization, 'b--', 'linewidth', 3);

% Label the obs likelihood
h_text_likelihood = text(2.5, 0.57, 'Likelihood y observed given \lambda', 'fontsize', 20);
h_text_posterior = text(1.4, 0.82, 'Posterior', 'fontsize', 20);
pause;
print -depsc s20f06.eps;

% Now zoom in on the bottom part
axis([1 2 0 2]);

for i = 1:num_points
   diff(i) = post_lam_val(i) / normalization - normpdf(lambda(i), lam_prior_mean, lam_prior_sd);
end

pause
print -depsc s20f07.eps;

% Turn off posterior label
set(h_text_likelihood, 'visible', 'off');
set(h_text_posterior, 'visible', 'off');

% Blow up the difference between prior and posterior for a minute
h_dif = plot(lambda, diff, 'm', 'linewidth', 3);
axis([1 2 -0.015 0.015]);
save_ticks = get(gca, 'ytick');
set(gca, 'ytick', [-0.01, 0, 0.01]);
h_diff_text = text(1.1, 0.01, '\lambda: Posterior - Prior', 'fontsize', 20);
h_diff_text2 = text(1.5, -0.01, 'Max density shifted to right', 'fontsize', 20);
grid on;

pause
print -depsc s20f08.eps;

% Reset the ticks
set(gca, 'ytick', save_ticks);
set(h_diff_text, 'visible', 'off');
set(h_diff_text2, 'visible', 'off');
axis([1 2 0 2]);
text(1.05, 1.7, 'Find Max by search', 'fontsize', 20);
text(1.35, 0.85, 'Max is new \lambda mean', 'fontsize', 20);


pause
print -depsc s20f09.eps;

