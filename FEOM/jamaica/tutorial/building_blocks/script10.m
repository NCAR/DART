% Finally, a kernel filter for some fun at maintaining distributions

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

axis([-4 4 0 0.85]);
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
set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');
print -depsc s10f01.eps;

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
product = prior .* obs;
h_prior_plot = plot(x, prior, 'g', 'linewidth', 3)
set(h_prior_plot, 'color', [0 0.73 0]);
h_prior_text = text(-3.0, 0.3, 'Prior PDF', 'fontsize', 24);
pause;
print -depsc s10f02.eps;

h_obs_plot = plot(x, obs, 'r', 'linewidth', 3)
h_obs_text = text(1.6, 0.4, 'Obs. Likelihood', 'fontsize', 24);
pause;
print -depsc s10f03.eps;

% Next, do a kernel fit around each ensemble member;
% Warning, kernel width is set to narrow here; not as in DART
set(h_prior_plot, 'linewidth', 1.0);
set(h_prior_text, 'visible', 'off');
set(p_text, 'visible', 'off');

% Label the kernels
h_kernel_text = text(-3.9, 0.72, 'Example Kernels: Half as Wide as Prior PDF', 'fontsize', 24);
for i = 1:5
   prior_kernel(i, :) = normpdf(x, x_prior(i), 1.2 / 2);
   h_prior_k(i) = plot(x, prior_kernel(i, :), 'g', 'linewidth', 3);
   set(h_prior_k(i), 'color', [0 0.73 0]);
   pause
   hhh = figure(1);
   outfile = ['s10f0', num2str(3 + i), '.eps'];
   print (hhh, '-depsc', outfile);
end

% Turn off the original fit
set(h_prior_plot, 'visible', 'off');
% Plot the normalized sum of the kernels (sum / 5 here) and lighten the kernels
for i = 1:5
   set(h_prior_k(i), 'linewidth', 1);
end
prior_sum = prior_kernel(1, :);
for i = 2:5
   prior_sum = prior_sum + prior_kernel(i, :);
end
prior_sum = prior_sum / 5;
h_prior_sum = plot(x, prior_sum, 'g', 'linewidth', 3);
set(h_prior_sum, 'color', [0 0.73 0]);
h_sum_text = text(-3.5, 0.43, 'Normalized Sum of Kernels', 'fontsize', 24);


pause
print -depsc s10f09.eps;

% Turn off the prior_sum for now
set(h_prior_sum, 'visible', 'off');
set(h_kernel_text, 'visible', 'off');
set(h_sum_text, 'visible', 'off');

pause
print -depsc s10f10.eps;

% Now do the products one by one; don't forget that there's now a weight associated
for i = 1:5
   % Highlight the appropriate prior kernel
   set(h_prior_k(i), 'linewidth', 3); 



   % Compute the shape of the product and plot
   product = prior_kernel(i, :) .* obs;
   pb = (sum(product) * 0.01); 
   posterior_kernel(i, :) = product./pb;
   h_pre_k(i) = plot(x, posterior_kernel(i, :), 'b', 'linewidth', 3)

   pause
   hhh = figure(1);
   outfile = ['s10f', num2str(9 + 2*i), '.eps'];
   print (hhh, '-depsc', outfile);

   % Turn down the unweighted product
   set(h_pre_k(i), 'visible', 'off');
   % Compute the weight for this kernel
   [temp, newmean_pos] = max(posterior_kernel(i, :));
   newmean = -5.0 + (newmean_pos - 1) * 0.01;
   newcov = 1. / (1. / (1.2 / 2) + 1. / 0.8);
   weight(i) = exp(-0.5 * (x_prior(i)^2 / (1.2 / 2) + 1^2 / 0.8 - newmean^2 / newcov));

   % Plot the weighted posterior kernel
   posterior_kernel(i, :) = posterior_kernel(i, :) * weight(i);
   h_posterior_k(i) = plot(x, posterior_kernel(i, :), 'b', 'linewidth', 3)

   pause 
   hhh = figure(1);
   outfile = ['s10f', num2str(10 + 2*i), '.eps'];
   print (hhh, '-depsc', outfile);
   
% Unhighlight previous prior and latest posterior
   set(h_prior_k(i), 'linewidth', 1); 
   set(h_posterior_k(i), 'linewidth', 1); 
end

% Turn off all the prior kernels
for i = 1:5
   set(h_prior_k(i), 'visible', 'off');
end

% Now plot the sum of the updated kernels and normalize
posterior_sum = posterior_kernel(1, :);
for i = 2:5
   posterior_sum = posterior_sum + posterior_kernel(i, :);
end
% Need to get integrated value of the observed value, y2                                 
pb = (sum(posterior_sum) * 0.01);                                                              
posterior_sum = posterior_sum ./ pb;               
h_posterior_sum = plot(x, posterior_sum, 'b', 'linewidth', 3);
h_post_sum_text = text(-3.9, 0.55, 'Normalized Sum of Posteriors', 'fontsize', 24);

pause
print -depsc s10f21.eps;

% Getting new sample is dicey, only easy way is random sampling 
% Would use weights for real; don't want to code this in matlab right now.
x_posterior = [1.02, -0.6 0.23 0.78 0.6];
y_posterior = [0.2 0.2 0.2 0.2 0.2];
h_f_plot = plot(x_posterior, y_posterior, 'b*')
set(h_f_plot, 'markersize', 18)
set(h_f_plot, 'linewidth', 3)
p_f_text = text(-2.0, 0.30, 'Posterior Ensemble', 'fontsize', 24);
print -depsc s10f22.eps;

