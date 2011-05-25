%% DART:script21 With ensemble methods, prior is available only as a 'random' sample

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

% Generate a random sample of the prior
x_prior = normrnd(-1.0, 1.2, 5, 1);
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
%text(-2.0, 0.1, 'Prior Ensemble', 'fontsize', 24);
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
% PRINTING IS SCREWED UP WITH CURRENT MATLAB VERSION, DO MANUALLY                        
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
hhh = plot(x, prior, 'g--', 'linewidth', 3)
set(hhh, 'color', [0 0.73 0]);
h_prior_text = text(-3.0, 0.3, 'Prior PDF', 'fontsize', 24);

plot(x, obs, 'r', 'linewidth', 3)
text(1.6, 0.4, 'Obs. Likelihood', 'fontsize', 24);
%h_text = text(-3.8, 0.56, '1. Compute updated inflation, lamda, distribution.', 'fontsize', 24);
pause;
print -depsc s21f01.eps;                       

set(h_prior_text, 'visible', 'off');

% Inflate prior by this amount
lamda = 1.5
x_prior_inf = sqrt(lamda) * (x_prior - mean(x_prior)) + mean(x_prior);
y_prior_inf = y_prior + 0.05;
h_plot = plot(x_prior_inf, y_prior_inf, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
set(h_plot, 'color', [0 0.73 0]);

prior_inf = normpdf(x, -1.0, 1.2 * sqrt(lamda));
product_inf = prior_inf .* obs;
h_inf = plot(x, prior_inf, 'g', 'linewidth', 3)
set(h_inf, 'color', [0 0.73 0]);

%set(h_text, 'visible', 'off');
%h_text = text(-3.8, 0.56, '2. Inflate ensemble with updated lamda mean.', 'fontsize', 24);

pause
print -depsc s21f02.eps;                       

% Need to get integrated value of the observed value, y2
pb = (sum(product_inf) * 0.01); 
posterior = product_inf./pb;
plot(x, posterior, 'b', 'linewidth', 3)
% Also get the updated ensemble (use EAF)
% First, need to compute the mean and variance of the blue
posterior_var = 1.0 / (1.0 / (lamda * 1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_mean = posterior_var * (-1.0 / (lamda * 1.2 * 1.2) + 1.0 / (0.8 * 0.8));
posterior_sd = sqrt(posterior_var);

% First shift the mean to have posterior mean
x_posterior = x_prior_inf - mean(x_prior_inf) + posterior_mean
y_posterior = y_prior + 0.15;

% Now adjust the variance by linear contraction
x_final = (x_posterior - posterior_mean) * sqrt(posterior_var / var(x_posterior)) + posterior_mean;
y_final = y_posterior + 0.1;
h_plot = plot(x_final, y_final, 'b*', 'markersize', 18, 'linewidth', 3)
%set(h_text, 'visible', 'off');
%h_text = text(-3.8, 0.56, '3. Compute posterior.', 'fontsize', 24);


pause
print -depsc s21f03.eps;                       

% Now show the total increments
[temp, indx] = sort(x_prior);
for i = 1:5
   xx = [x_prior(indx(i)), x_final(indx(i))];
   yy = [0.14 + 0.02*i,  0.14 + 0.02 * i];
   plot (xx, yy, 'b', 'linewidth', 3);
   % Draw a thin vertical to the original prior
   xx = [x_prior(indx(i)), x_prior(indx(i))];
   yy = [0.14 + 0.02*i, 0.05];
   plot(xx, yy, 'b');
   % Draw a thin vertical to the posterior
   xx = [x_final(indx(i)), x_final(indx(i))];
   yy = [0.14 + 0.02*i, y_final(i)];
   plot(xx, yy, 'b');
end

%set(h_text, 'visible', 'off');
%h_text = text(-3.8, 0.56, '4. Get increments from original ensemble.', 'fontsize', 24);

print -depsc s21f04.eps

