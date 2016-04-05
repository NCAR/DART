%% DART:script11

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

% Printing gums up things here, have to use print_count to get later figures
print_count = 13;

% set this to 0 if you do not want to pause before each figure
wait = 1;

% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 6.5 6.5]);
% Set the shape of the plot box
%pbaspect([2.4 1 1]);

% Take a look at using multiple axes here
r1 = subplot(2, 2, 1);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
h_y_label = ylabel('Unobserved State Variable');
r2 = subplot(2, 2, 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
r3 = subplot(2, 2, 4);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
xlabel('Observed Variable');


% Position the marginal and joint plot boxes
set(r1, 'Position', [0.14, 0.23, 0.10, 0.6550]);
set(r2, 'Position', [0.25, 0.23, 0.6550, 0.6550]);
set(r3, 'Position', [0.25 0.1100 0.6550 0.1000]);

% Turn off the unwanted tick labels
set(r2, 'Xticklabel', []);
set(r2, 'yticklabel', []);
set(r1, 'Xticklabel', []);
set(r3, 'yticklabel', []);

subplot(r2);
% First look at impacts on a second 'unobserved' variable
% Need to be able to generate a two dimensional random gaussian draw
ens_size = 5;
for i = 1:ens_size
   means = [1 4];
   correlation = 0.8;
   variance = [4.0 0.5];
   covariance = correlation * sqrt(variance(1)) * sqrt(variance(2));
   % Use Knuth method to generate
   a11 = sqrt(variance(1));
   a21 = covariance / a11;
   a22 = sqrt(variance(2) - a21^2);
   % Two base independent gaussian deviates
   randx1 = normrnd(0, 1);
   randx2 = normrnd(0, 1);
   % Use these to generate correlated
   rnum(1, i) = means(1) + a11 * randx1;
   rnum(2, i) = means(2) + a21 * randx1 + a22 * randx2;
end


r_joint_prior = plot(rnum(1, :), rnum(2, :), 'g*');
set(r_joint_prior, 'markersize', 12);
set(r_joint_prior, 'linewidth', 2);
set(r_joint_prior, 'color', [0 0.73 0]);
grid on;
set(r2, 'Xticklabel', []);
set(r2, 'yticklabel', []);
hold on;

% Compute the regression line,
s_correl = corrcoef(rnum(1, :), rnum(2, :));
s_var(1) = var(rnum(1, :));
s_var(2) = var(rnum(2, :));
s_covar = s_correl * sqrt(s_var(1)) * sqrt(s_var(2));
reg_coef = s_covar(1, 2) / s_var(1);



% Plot the marginals in the wing plots
subplot(r3);
yu(1:ens_size) = 0.1;
h_prior_marg = plot(rnum(1, :), yu, '*g');
set(h_prior_marg, 'markersize', 12)
set(h_prior_marg, 'linewidth', 2);
set(h_prior_marg, 'color', [0 0.73 0]);
% Use the x axis limits from the main plot
axis([get(r2, 'xlim')  [0 1]]);
set(r3, 'yticklabel', []);
xlabel('Observed Variable');

% Plot the marginals in the wing plots
subplot(r1);
xu(1:ens_size) = 0.1;
h_state_marg = plot(xu, rnum(2, :), '*g');
set(h_state_marg, 'color', [0 0.73 0]);
set(h_state_marg, 'markersize', 12)
set(h_state_marg, 'linewidth', 2);
axis([[0 1] get(r2, 'ylim')]);
set(r1, 'xticklabel', []);
h_y_label = ylabel('Unobserved State Variable');

if (wait > 0), pause, end;

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s11f01.eps;

set(r1, 'Position', [0.14, 0.73, 0.10, 0.1550]);
set(r2, 'Position', [0.25, 0.73, 0.6550, 0.1550]);
set(r3, 'Position', [0.25 0.1100 0.6550 0.6000]);

h_y_label = ylabel('Unobs.');
%set(h_y_label, 'visible', 'off');


subplot(r3);
hold on;
% Do a simple update scheme in the marginal window
% Generate a random sample of the prior
x_prior = rnum(1, :);

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, mean(x_prior), std(x_prior));
obs = normpdf(x, 2.3, 0.6);
product = prior .* obs;
h_prior = plot(x, prior, 'g', 'linewidth', 3)
set(h_prior, 'color', [0 0.73 0]);
%text(-3.0, 0.3, 'Prior PDF', 'fontsize', 24);
if (wait > 0), pause, end;
print -depsc s11f02.eps;

h_obs = plot(x, obs, 'r', 'linewidth', 3)
%text(1.6, 0.4, 'Obs. Likelihood', 'fontsize', 24);
if (wait > 0), pause, end;
print -depsc s11f03.eps;

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01);
posterior = product./pb;
h_posterior = plot(x, posterior, 'b', 'linewidth', 3)
%text(-1.9, 0.55, 'Posterior PDF', 'fontsize', 24);

if (wait > 0), pause, end;
print -depsc s11f04.eps;

% Plot a random sample from the update (cheat and take it from
% obs which will be very similar here).
x_posterior = normrnd(2.3, 0.6, 5, 1);
y_posterior = [0.9 0.9 0.9 0.9 0.9];
h_post_plot = plot(x_posterior, y_posterior, 'b*');
set(h_post_plot, 'markersize', 12)
set(h_post_plot, 'linewidth', 2);

if (wait > 0), pause, end;
print -depsc s11f05.eps;

% Plot the increment vectors
h_inc = text(-1.5, 0.4, 'Increments', 'fontsize', 24);
for i = 1:5
   arr1 = [x_prior(i), 0.15 + 0.175 * (i-1)];
   arr2 = [x_posterior(i), 0.15 + 0.175 * (i-1)];
   h_arrow(i) = arrow(arr1, arr2, 4.0);
   set(h_arrow(i), 'linewidth', 3);
   if (wait > 0), pause, end;
   outfile = ['s11f0', num2str(5 + i), '.eps'];
   print (gcf, '-depsc', outfile);

end

% Turn off the update continuous distributions
set(h_prior, 'visible', 'off');
set(h_posterior, 'visible', 'off');
set(h_obs, 'visible', 'off');

%Turn off the old prior and posterior asterisks and redraw them at the arrows tips
set(h_post_plot, 'visible', 'off');
set(h_prior_marg, 'visible', 'off');
hold on;
for i = 1:5
   h_prior(i) = plot(x_prior(i), 0.15 + 0.175 * (i -1), 'g*');
   set(h_prior(i), 'markersize', 12)
   set(h_prior(i), 'linewidth', 2);
   set(h_prior(i), 'color', [0 0.73 0]);
   h_posterior(i) = plot(x_posterior(i), 0.15 + 0.175 * (i - 1), 'b*');
   set(h_posterior(i), 'markersize', 12)
   set(h_posterior(i), 'linewidth', 2);
end

if (wait > 0), pause, end;
print -depsc s11f11.eps;

% Position the marginal and joint plot boxes
set(r1, 'Position', [0.14, 0.23, 0.10, 0.6550]);
set(r2, 'Position', [0.25, 0.23, 0.6550, 0.6550]);
set(r3, 'Position', [0.25 0.1100 0.6550 0.1000]);

subplot(r1);
h_y_label = ylabel('Unobserved State Variable');
%set(h_y_label, 'visible', 'on');

%Plot point at mean of distibution in black
s_mean(1) = mean(rnum(1, :));
s_mean(2) = mean(rnum(2, :));
%plot(s_mean(1), s_mean(2), 'r*');
%Move the label for the squeeze
set(h_inc, 'position', [-1.5 0.82 0]);

if (wait > 0), pause, end;
print -depsc s11f12.eps;

% Plot best fit line through this mess; just uses the standard deviations?
subplot(r2);
x_point(1) = min(rnum(1, :));
x_point(2) = max(rnum(1, :));
slope = sqrt(s_var(1)) * sqrt(s_var(2)) / s_var(1);
y_point(1) = s_mean(2) + (x_point(1) - s_mean(1)) * slope;
y_point(2) = s_mean(2) + (x_point(2) - s_mean(1)) * slope;
plot(x_point, y_point, 'r', 'linewidth', 3);

if (wait > 0), pause, end;
if(print_count == 13)
   print -depsc s11f13.eps;
   print_count = print_count+1;
end

% Plot the un-normalized increment vectors
% Turn all the incrment lines down
subplot(r3)
for i = 1:5
   set(h_arrow(i), 'linewidth', 1);
end
for i = 1:5
   subplot(r3)
   set(h_arrow(i), 'linewidth', 3);
   subplot(r2)
   x_inc = [x_prior(i), x_posterior(i)];
   y_delta(i) = (x_posterior(i) - x_prior(i)) * slope
   y_inc = [rnum(2, i), rnum(2, i) + (x_posterior(i) - x_prior(i)) * slope];
   h_slope(i) = plot(x_inc, y_inc, 'b', 'linewidth', 3);
   if (wait > 0), pause, end;
   outfile = ['s11f', num2str(13 + i), '.eps'];
   if(print_count == 13 + i)
      print (gcf, '-depsc', outfile);
      print_count = print_count+1;
   end

   subplot(r3)
   set(h_arrow(i), 'linewidth', 1);
   subplot(r2)
   set(h_slope(i), 'linewidth', 1);
end

% Next, do the projection onto the state variable
% First, move the prior state variables to have varied positions
set(h_state_marg, 'visible', 'off');
subplot(r1);
hold on;
for i = 1:5
   h_state_prior(i) = plot(0.15 + 0.175 * (i - 1), rnum(2, i), 'g*');
   set(h_state_prior(i), 'markersize', 12)
   set(h_state_prior(i), 'linewidth', 2);
   set(h_state_prior(i), 'color', [0 0.73 0]);
end

for i = 1:5
   subplot(r3)
   set(h_arrow(i), 'linewidth', 3);
   subplot(r2)
   set(h_slope(i), 'linewidth', 3);
   % First, plot the straight projection in light blue
   xv = [0.15 + 0.175 * (i - 1), 0.15 + 0.175 * (i - 1)];
   yv = [rnum(2, i), rnum(2, i) + y_delta(i)];
   subplot(r1);
   h_full_dist(i) = plot(xv, yv);
   yvc = [rnum(2, i), rnum(2, i) + y_delta(i) * s_correl(1, 2)];
   h_part_dist(i) = plot(xv, yvc, 'linewidth', 3);
   if (wait > 0), pause, end;
   outfile = ['s11f', num2str(18 + i), '.eps'];
   if(print_count == 18 + i)
      print (gcf, '-depsc', outfile);
      print_count = print_count+1;
   end

   subplot(r3)
   set(h_arrow(i), 'linewidth', 1);
   subplot(r2)
   set(h_slope(i), 'linewidth', 1);
end

% Turn off the full distance lines and put in updated asterisksA
for i = 1:5
   subplot(r1);
   set(h_full_dist(i), 'visible', 'off');
   final_state(i) = rnum(2, i) + y_delta(i) * s_correl(1, 2);
   h_final_state(i) = plot(0.15 + 0.175 * (i - 1), final_state(i), 'b*');
   set(h_final_state(i), 'markersize', 12)
   set(h_final_state(i), 'linewidth', 2);
   subplot(r3)
   set(h_arrow(i), 'linewidth', 3);
   subplot(r2)
   set(h_slope(i), 'linewidth', 3);
end

% Now take a quick peak at what's going on in the state variable
set(r1, 'Position', [0.14, 0.23, 0.6000, 0.6550]);
set(r2, 'Position', [0.75, 0.23, 0.1550, 0.6550]);
set(r3, 'Position', [0.75 0.1100 0.1550 0.1000]);

% Move the Green asterisks and blue asterisks to one line, then do fit
set(h_inc, 'visible', 'off');
subplot(r3);
xlabel('Obs.');
subplot(r1);
for i = 1:5
   set(h_part_dist(i), 'visible', 'off');
   set(h_final_state(i), 'xdata', [0.15]);
   set(h_state_prior(i), 'xdata', [0.05]);
end

if (wait > 0), pause, end;
if(print_count == 24)
   print -depsc s11f24.eps;
   print_count = print_count+1;
end

% Get 1000 points from bottom of axis to top of axis
y_limits = get(gca, 'ylim');
y_final_pts = y_limits(1) : 0.01: y_limits(2);
% Plot the prior fit
prior_mean = mean(rnum(2, :));
prior_std = std(rnum(2, :));
x_final_pts = normpdf(y_final_pts, prior_mean, prior_std);
hhh = plot(x_final_pts, y_final_pts, 'g', 'linewidth', 3);
set(hhh, 'color', [0 0.73 0]);
text(0.4, 3.3, 'Prior State Fit', 'fontsize', 24);


if (wait > 0), pause, end;
if(print_count == 25)
   print -depsc s11f25.eps;
   print_count = print_count+1;
end

% Plot the posterior fit
final_mean = mean(final_state);
final_std = std(final_state);
x_final_pts = normpdf(y_final_pts, final_mean, final_std);
plot(x_final_pts, y_final_pts, 'linewidth', 3);
text(0.55, 4.9, 'Posterior Fit', 'fontsize', 24);
if(print_count == 26)
   print -depsc s11f26.eps;
end

