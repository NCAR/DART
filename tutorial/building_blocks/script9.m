%% DART:script9 

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
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
randn('state', 6);
% Classical perturbed obs ensemble Kalman filter

x_prior = [-2.2 -1.8 -1.4 1.8 2.2];
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior = (x_prior - mean(x_prior)) * 1.2 / std(x_prior) - 1.0

% Plot the prior distribution
y_prior = [0.05, 0.05, 0.05, 0.05, 0.05];
for i = 1:5
   h_prior_plot(i) = plot(x_prior(i), y_prior(i), 'g*');
   hold on;
   set(h_prior_plot(i), 'markersize', 18);
   set(h_prior_plot(i), 'linewidth', 3);
   set(h_prior_plot(i), 'color', [0 0.73 0]);
end

axis([-4 4 0 0.6]);
set(gca, 'fontsize', 24);
p_text = text(-3.8, 0.1, 'Prior Ensemble', 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');
pause

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s09f01.eps;

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
product = prior .* obs;
h_prior = plot(x, prior, 'g', 'linewidth', 3)
set(h_prior, 'color', [0 0.73 0]);
%h_prior_text = text(-3.0, 0.3, 'Prior PDF', 'fontsize', 24);
pause;
print -depsc s09f02.eps;

h_obs_base = plot(x, obs, 'r', 'linewidth', 3)
h_obs_text = text(1.6, 0.4, 'Obs. Likelihood', 'fontsize', 24);
pause;
print -depsc s09f03.eps;

% Pick a random draw from the obs. likelihood 
x_obs = normrnd(1.0, 0.8, 5, 1);
y_obs = [0.20 0.2 0.2 0.2 0.2];

% Draw a line from each prior member to random obs
rd_text = text(0.2, 0.25, 'Random Draws from Obs.', 'fontsize', 24);
for i = 1:5
   h_obs_plot(i) = plot(x_obs(i), y_obs(i), 'r*')
   set(h_obs_plot(i), 'markersize', 18)
   set(h_obs_plot(i), 'linewidth', 3)
   xu = [x_prior(i), x_obs(i)];
   yu = [y_prior(i), y_obs(i)];
   h_obs_line(i) = plot(xu, yu, 'k')

   pause;
   hhh = figure(1);
   outfile = ['s09f0', num2str(3 + i), '.eps'];
   print (hhh, '-depsc', outfile);

end

% Make all but the first pair disappear
for i = 2:5
   set(h_obs_plot(i), 'visible', 'off');
   set(h_obs_line(i), 'visible', 'off');
   set(h_prior_plot(i), 'visible', 'off');
end
set(h_obs_text, 'visible', 'off');
set(rd_text, 'visible', 'off');
set(p_text, 'visible', 'off');

pause;
print -depsc s09f09.eps;

% Shift a copy of the prior variance over the first prior
prior_shift = normpdf(x, x_prior(1), 1.2);
h_prior_shift = plot(x, prior_shift, 'g', 'linewidth', 3)
set(h_prior_shift, 'color', [0 0.73 0]);
% Draw a dashed line marking the mean to the point
xm = [x_prior(1) x_prior(1)];
ym = [y_prior(1) max(prior)];
h_prior_vert = plot(xm, ym, 'g--', 'linewidth', 2);
set(h_prior_vert, 'color', [0 0.73 0]);
% Lighten up the prior pdf original and dot it
set(h_prior, 'linewidth', 1);
set(h_prior, 'linestyle', '--');
pause
print -depsc s09f10.eps;

% Shift a copy of obs variance over the first perturbed obs
obs_shift = normpdf(x, x_obs(1), 0.8);
h_obs_shift = plot(x, obs_shift, 'r', 'linewidth', 3);
% Draw a dashed line marking the mean to the point
xm = [x_obs(1) x_obs(1)];
ym = [y_obs(1) max(obs)];
h_obs_vert = plot(xm, ym, 'r--', 'linewidth', 2);
% Lighten up the obs pdf original
set(h_obs_base, 'linewidth', 1);
set(h_obs_base, 'linestyle', '--');
pause;
print -depsc s09f11.eps;


% Need to get integrated value of the product of shifted
product = prior_shift .* obs_shift;
pb = (sum(product) * 0.01); 
posterior = product./pb;
h_posterior_shift = plot(x, posterior, 'b', 'linewidth', 3)
h_post_text = text(-2.2, 0.55, 'Posterior PDF', 'fontsize', 24);
% Also plot a point at the mean value
[post_max_y post_max_index] = max(posterior);
post_max_x = -5.0 + (post_max_index - 1) * 0.01;
posterior_x_ens(1) = post_max_x;
h_post_vert = plot([post_max_x post_max_x], [0.3, post_max_y], 'b--', 'linewidth', 2);
h_post_shift_mean = plot(post_max_x, 0.3, 'b*');
set(h_post_shift_mean, 'markersize', 18)
set(h_post_shift_mean, 'linewidth', 3)
pause;
print -depsc s09f12.eps;

% Make the shifted distributions and the mean vertical lines disappear
set(h_obs_shift, 'visible', 'off');
set(h_prior_shift, 'visible', 'off');
set(h_posterior_shift, 'visible', 'off');
set(h_post_vert, 'visible', 'off');
set(h_prior_vert, 'visible', 'off');
set(h_obs_vert, 'visible', 'off');
set(h_obs_line(1), 'visible', 'off');
set(h_post_text, 'visible', 'off');

pause
print -depsc s09f13.eps;

% Loop through the remaining prior-obs pairs and get updated
for i = 2:5
   % Plot the prior and obs pair and turn on the line between them
   set(h_obs_line(i), 'visible', 'on');
   set(h_obs_plot(i), 'visible', 'on');
   set(h_prior_plot(i), 'visible', 'on');

   % Plot the shifted distributions but be prepared to turn them back off
   % Shift a copy of the prior variance over the first prior
   prior_shift = normpdf(x, x_prior(i), 1.2);
   h_prior_shift = plot(x, prior_shift, 'g', 'linewidth', 3)
   set(h_prior_shift, 'color', [0 0.73 0]);
   % Draw a dashed line marking the mean to the point
   xm = [x_prior(i) x_prior(i)];
   ym = [y_prior(i) max(prior)];
   h_prior_vert = plot(xm, ym, 'g--', 'linewidth', 2);
    set(h_prior_vert, 'color', [0 0.73 0]);

   % Shift a copy of obs variance over the first perturbed obs
   obs_shift = normpdf(x, x_obs(i), 0.8);
   h_obs_shift = plot(x, obs_shift, 'r', 'linewidth', 3);
   % Draw a dashed line marking the mean to the point
   xm = [x_obs(i) x_obs(i)];
   ym = [y_obs(i) max(obs)];
   h_obs_vert = plot(xm, ym, 'r--', 'linewidth', 2);

   % Need to get integrated value of the product of shifted
   product = prior_shift .* obs_shift;
   pb = (sum(product) * 0.01); 
   posterior = product./pb;
   h_posterior_shift = plot(x, posterior, 'b', 'linewidth', 3)
   % Also plot a point at the mean value
   [post_max_y post_max_index] = max(posterior);
   post_max_x = -5.0 + (post_max_index - 1) * 0.01;
   posterior_x_ens(i) = post_max_x;
   h_post_vert = plot([post_max_x post_max_x], [0.3, post_max_y], 'b--', 'linewidth', 2);
   h_post_shift_mean = plot(post_max_x, 0.3, 'b*');
   set(h_post_shift_mean, 'markersize', 18)
   set(h_post_shift_mean, 'linewidth', 3)

   pause
   hhh = figure(1);
   outfile = ['s09f', num2str(12 + i), '.eps'];
   print (hhh, '-depsc', outfile);

   % Turn off the line in between this pair
   set(h_obs_line(i), 'visible', 'off');
   % Turn off the shifted distributions
   set(h_prior_shift, 'visible', 'off');
   set(h_prior_vert, 'visible', 'off')
   set(h_obs_shift, 'visible', 'off');
   set(h_obs_vert, 'visible', 'off')
   set(h_posterior_shift, 'visible', 'off');
   set(h_post_vert, 'visible', 'off')
end

% Finally, turn on lines between the initial greens and the posteriors
for i = 1:5
   xu = [x_prior(i) posterior_x_ens(i)];
   yu = [0.05 0.3];
   plot(xu, yu, 'k');
end

print -depsc s09f18.eps;

