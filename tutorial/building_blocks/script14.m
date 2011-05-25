%% DART:script14 Demonstration of how an unrelated observation evolves a State variable due to sampling error in the regression

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
set(gcf, 'position', [1 1 6.5 6.5]);
% Set the shape of the plot box
%pbaspect([2.4 1 1]);

% Take a look at using multiple axes here
r1 = subplot(2, 2, 1);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
ylabel('Unobserved State Variable');
r2 = subplot(2, 2, 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
r3 = subplot(2, 2, 4);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
xlabel('Observed Variable');


% Position the marginal and joint plot boxes
set(r1, 'Position', [0.14, 0.23, 0.25, 0.6550]);
set(r2, 'Position', [0.40, 0.23, 0.5050, 0.6550]);
set(r3, 'Position', [0.40 0.1100 0.5050 0.1000]);

% Turn off the unwanted tick labels
set(r2, 'Xticklabel', []);
set(r2, 'yticklabel', []);
set(r1, 'Xticklabel', []);
set(r3, 'yticklabel', []);

subplot(r2);
% Plot something invisible to make first panel look consistent
h_temp = plot([0 0], [1 1]);
set(gca, 'linewidth', 2);
axis([-3 3 -3 3]);
set(h_temp, 'visible', 'off');
set(r2, 'Xticklabel', []);
set(r2, 'yticklabel', []);
hold on;

subplot(r3);
% Plot something invisible to make first panel look consistent
h_temp = plot([0 0], [1 1]);
set(gca, 'linewidth', 2);
axis([ -3 3 0 1]);
set(h_temp, 'visible', 'off');
set(r3, 'yticklabel', []);
xlabel('Observed Variable');
hold on;


% First look at impacts on a second 'unobserved' variable
% Need to be able to generate a two dimensional random gaussian draw
ens_size = 20;
state = normrnd(0, 1, ens_size, 1);
subplot(r1);
state_x(1:ens_size) = 0.05;
h_state_prior = plot(state_x, state, 'g*');
set(h_state_prior, 'color', [0 0.73 0]);
hold on;
ylabel('Unobserved State Variable');
set(gca, 'linewidth', 2);
sd_label = ['SD=', sprintf('%0.2f', std(state))];
text(0.15, 2.80, sd_label, 'fontsize', 20, 'color', [0 0.73 0]);
mean_label = ['MN=', sprintf('%0.2f', mean(state))];
text(0.15, 2.45, mean_label, 'fontsize', 20, 'color', [0 0.73 0]);
set(h_state_prior, 'markersize', 12);
set(h_state_prior, 'linewidth', 2)
axis([0 1 -3 3]);
set(r1, 'Xticklabel', []);

pause
                                                                                         
% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s14f01.eps;

% Observe an unrelated variable repeatedly
for i = 1:121
   obs = normrnd(0, 1, ens_size, 1);

   % Observe a random draw
   obs_val = normrnd(0, 1);
   % Compute prior obs variance
   prior_obs_var = var(obs);
   prior_obs_mean = mean(obs);
   % Compute updated variance
   update_var = 1.0 / (1.0 / 1.0 + 1.0 / prior_obs_var);
   % Updated mean is
   new_mean = update_var * (prior_obs_mean / prior_obs_var + obs_val / 1.0);
   posterior_obs = (obs - mean(obs)) * sqrt(update_var / 1.0) + new_mean;
   obs_inc = posterior_obs - obs;

   % Compute the regression coefficient
   s_correl = corrcoef(state, obs);
   s_var(1) = var(obs);
   s_var(2) = var(state);
   s_covar = s_correl(1, 2) * sqrt(s_var(1)) * sqrt(s_var(2));
   reg_coef = s_covar / s_var(1);

   % Compute posterior state
   post_state = state + reg_coef * obs_inc;
   output_interval = 20;
   if mod(i, output_interval) == 1
      % Plot the posterior observation
      subplot(r3); 
      h_obs_prior = plot(obs, state_x, 'g*');
      set(h_obs_prior, 'markersize', 12);
      set(h_obs_prior, 'linewidth', 2)
      set(h_obs_prior, 'color', [0 0.73 0]);
      xlabel('Observed Variable');
      set(r3, 'yticklabel', []);
      hold on;
      set(gca, 'linewidth', 2);
      axis([ -3 3 0 1]);
      % Plot the observation likelihood
      xlike = -3:0.01:3;
      ylike = normpdf(xlike, obs_val, 1.0);
      h_obs_like = plot(xlike, 2.3*ylike, 'r', 'linewidth', 2);
      % Plot the joint prior distribution
      subplot(r2);
      h_joint_prior = plot(obs, state, 'g*');
      set(h_joint_prior, 'color', [0 0.73 0]);
      hold on;
      set(gca, 'linewidth', 2);
      set(r2, 'Xticklabel', []);
      set(r2, 'yticklabel', []);
      axis([-3 3 -3 3]);
      set(h_joint_prior, 'markersize', 12);
      set(h_joint_prior, 'linewidth', 2)
      % Plot the observation variable posterior
      subplot(r3);
      hold on;
      h_obs_post = plot(posterior_obs, state_x + 0.3, 'b*');
      set(h_obs_post, 'markersize', 12);
      set(h_obs_post, 'linewidth', 2)
      
      % Plot the joint posterior distribution
      subplot(r2);
      h_joint_post = plot(posterior_obs, state, 'b*');
      set(h_joint_post, 'markersize', 12);
      set(h_joint_post, 'linewidth', 2)
      % Label the step number and the correlation
      step_label = ['After Obs. ', num2str(i)];
      h_step = text(-2.8, 2.8, step_label, 'fontsize', 24);
      correl_label = ['Sample Correl. = ', sprintf('%0.2f', s_correl(1, 2))];
      h_correl = text( -2.8, -2.8, correl_label, 'fontsize', 20);

      % Plot the state variable update
      subplot(r1);
      hold on;
      h_state_post = plot(0.15 + state_x + 0.15 * i / output_interval, post_state, 'b*');
      set(h_state_post, 'markersize', 12);
      set(h_state_post, 'linewidth', 2);
      sd_label = ['SD=', sprintf('%0.2f', std(post_state))];
      h_sd = text(0.15, -2.50, sd_label, 'fontsize', 20, 'color', [0 0 1]);
      mean_label = ['MN=', sprintf('%0.2f', mean(post_state))];
      h_mean = text(0.15, -2.85, mean_label, 'fontsize', 20, 'color', [0 0 1]);
      pause;

      outname = ['s14f0', num2str(2 + (i - 1) / 20)]
      print(gcf, '-depsc', outname);

      % Clear out all the junk
      set(h_obs_post, 'visible', 'off');
      set(h_obs_prior, 'visible', 'off');
      set(h_joint_prior, 'visible', 'off');
      set(h_joint_post, 'visible', 'off');
      set(h_obs_like, 'visible', 'off');
      set(h_step, 'visible', 'off');
      set(h_correl, 'visible', 'off');
      set(h_sd, 'visible', 'off');
      set(h_mean, 'visible', 'off');
      
   end
   state = post_state;

end

