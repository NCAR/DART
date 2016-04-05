% Looking at effects of errors in least squares approximation

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

% Set up a non-linear relation between priors and proceed

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
ens_size = 5;
% Use the rnum array to hold something decidely not random
rnum(1:2, 1:5) = 0;
rnum(1, :) = -1.5: 1.25: 3.5;
rnum(2, :) = rnum(1, :) .^2 + 2;




r_joint_prior = plot(rnum(1, :), rnum(2, :), 'g*');
set(r_joint_prior, 'markersize', 12);
set(r_joint_prior, 'linewidth', 2);
set(r_joint_prior, 'color', [0 0.73 0]);
grid on;
set(r2, 'Xticklabel', []);
set(r2, 'yticklabel', []);
hold on;

% Plot the y = x^2 curve on which the relation actually lies
x_limits = get(gca, 'xlim');
x_pts = x_limits(1) : 0.01: x_limits(2);
y_pts = x_pts .^ 2 + 2;
plot(x_pts, y_pts, 'r', 'linewidth', 2);

pause
% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s12f01.eps;

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
set(h_prior_marg, 'color', [0 0.73 0]);
set(h_prior_marg, 'markersize', 12)
set(h_prior_marg, 'linewidth', 2);
% Use the x axis limits from the main plot
axis([get(r2, 'xlim')  [0 1]]);
set(r3, 'yticklabel', []);
xlabel('Observed Variable');

% Plot the marginals in the wing plots
subplot(r1);
xu(1:ens_size) = 0.1;
h_state_marg = plot(xu, rnum(2, :), '*g');
set(h_state_marg, 'markersize', 12)
set(h_state_marg, 'linewidth', 2);
set(h_state_marg, 'color', [0 0.73 0]);
axis([[0 1] get(r2, 'ylim')]);
set(r1, 'xticklabel', []);
h_y_label = ylabel('Unobserved State Variable');

pause;
print -depsc s12f02.eps;

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
%pause;                                                                                   
                                                                                         
h_obs = plot(x, obs, 'r', 'linewidth', 3)                                                        
% Need to get integrated value of the observed value, y2                                 
pb = (sum(product) * 0.01);                                                              
posterior = product./pb;                                                                 
h_posterior = plot(x, posterior, 'b', 'linewidth', 3)                                                  

% Plot a random sample from the update (cheat and take it from
% obs which will be very similar here).
x_posterior = normrnd(2.3, 0.6, 5, 1);
y_posterior = [0.9 0.9 0.9 0.9 0.9];
h_post_plot = plot(x_posterior, y_posterior, 'b*');
set(h_post_plot, 'markersize', 12)
set(h_post_plot, 'linewidth', 2);

pause
print -depsc s12f03.eps;

% Plot the increment vectors
h_inc = text(-1.5, 0.4, 'Increments', 'fontsize', 24);
for i = 1:5
   arr1 = [x_prior(i), 0.15 + 0.175 * (i-1)];
   arr2 = [x_posterior(i), 0.15 + 0.175 * (i-1)];
   h_arrow(i) = arrow(arr1, arr2, 4.0);
   set(h_arrow(i), 'linewidth', 3);
end

pause
print -depsc s12f04.eps;

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

% Position the marginal and joint plot boxes
set(r1, 'Position', [0.14, 0.23, 0.10, 0.6550]);
set(r2, 'Position', [0.25, 0.23, 0.6550, 0.6550]);
set(r3, 'Position', [0.25 0.1100 0.6550 0.1000]);
subplot(r1);
h_y_label = ylabel('Unobserved State Variable');
%set(h_y_label, 'visible', 'on');

subplot(r2);
%Plot point at mean of distibution in black
s_mean(1) = mean(rnum(1, :));
s_mean(2) = mean(rnum(2, :));
%plot(s_mean(1), s_mean(2), 'r*');
%Move the label for the squeeze
set(h_inc, 'position', [-1.5 0.82 0]);

pause
print -depsc s12f05.eps;

% Plot best fit line through this mess; just uses the standard deviations?
subplot(r2);
%x_point(1) = min(rnum(1, :));
x_point(1) = -1.4;
%x_point(2) = max(rnum(1, :));
x_point(2) = 4;
slope = sqrt(s_var(1)) * sqrt(s_var(2)) / s_var(1);
y_point(1) = s_mean(2) + (x_point(1) - s_mean(1)) * slope;
y_point(2) = s_mean(2) + (x_point(2) - s_mean(1)) * slope;
plot(x_point, y_point, 'r', 'linewidth', 3);

pause
print -depsc s12f06.eps;

% Plot the un-normalized increment vectors
% Turn all the incrment lines down
for i = 1:5
   set(h_arrow(i), 'linewidth', 1);
end
for i = 1:5
   set(h_arrow(i), 'linewidth', 3);
   x_inc = [x_prior(i), x_posterior(i)];
   y_delta(i) = (x_posterior(i) - x_prior(i)) * slope
   y_inc = [rnum(2, i), rnum(2, i) + (x_posterior(i) - x_prior(i)) * slope];
   h_slope(i) = plot(x_inc, y_inc, 'b', 'linewidth', 3);
end

pause
print -depsc s12f07.eps;

% Try sorting the increments instead; begin by sorting the posteriors
% Priors are already sorted in this instance so don't need to sort them
x_posterior = sort(x_posterior);

% Blow up the observed variable marginal box again
set(r1, 'Position', [0.14, 0.73, 0.10, 0.1550]);
set(r2, 'Position', [0.25, 0.73, 0.6550, 0.1550]);
set(r3, 'Position', [0.25 0.1100 0.6550 0.6000]);
subplot(r1);
h_y_label = ylabel('Unobs.');
%set(h_y_label, 'visible', 'off');

pause
print -depsc s12f08.eps;

% Replace with the sorted posteriors one by one
for i = 1:5
   set(h_posterior(i), 'xdata', x_posterior(i));
end

pause
print -depsc s12f09.eps;

subplot(r3);
% Now update the arrows
for i = 1:5
   set(h_arrow(i), 'visible', 'off');
   arr1 = [x_prior(i), 0.15 + 0.175 * (i-1)];
   arr2 = [x_posterior(i), 0.15 + 0.175 * (i-1)];
   h_arrow2(i) = arrow(arr1, arr2, 4.0);
   set(h_arrow2(i), 'linewidth', 3);
end

pause
print -depsc s12f10.eps;

% Back to base positioning
set(r1, 'Position', [0.14, 0.23, 0.10, 0.6550]);
set(r2, 'Position', [0.25, 0.23, 0.6550, 0.6550]);
set(r3, 'Position', [0.25 0.1100 0.6550 0.1000]);
subplot(r1);
h_y_label = ylabel('Unobserved State Variable');
%set(h_y_label, 'visible', 'on');

%Replace the old update vectors with the new ones one by one;
subplot(r2);
for i = 1:5
   % Put a star on to mark the new position
   y_delta(i) = (x_posterior(i) - x_prior(i)) * slope
   h_post2(i) = plot(x_posterior(i), rnum(2, i) + y_delta(i), 'b*');  
   set(h_post2(i), 'markersize', 12)
   set(h_post2(i), 'linewidth', 2);
end

pause
print -depsc s12f11.eps;

% Finally, just leave the new update vectors
for i = 1:5
   set(h_slope(i), 'visible', 'off');
   x_inc = [x_prior(i), x_posterior(i)];
   y_delta(i) = (x_posterior(i) - x_prior(i)) * slope
   y_inc = [rnum(2, i), rnum(2, i) + (x_posterior(i) - x_prior(i)) * slope];
   h_slope(i) = plot(x_inc, y_inc, 'b', 'linewidth', 3);
end
print -depsc s12f12.eps;

