%% DART:old4

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

%% Plots for a scalar update
figure(2)
x = -5:0.01:5;
prior = normpdf(x, -1.0);
obs = normpdf(x, 1.0);
product = prior .* obs;
plot(x, prior, 'g')
hold on;
plot(x, obs, 'r')
% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)



figure(4)
% Make a random sample with mean -1 and variance 1
prior_x = normrnd(-1, 1, 4, 1)
% Adjust them to have mean -1 and variance 1
prior_x = prior_x - mean(prior_x) -1;
sd = sqrt(var(prior_x))
prior_x = (prior_x - mean(prior_x)) * (1.0 / sd) + mean(prior_x)
prior_y = [0.05 0.05 0.05 0.05];
h = plot(prior_x, prior_y, 'g*') 
%set(h, 'fontsize', 24)
hold on;
x = -5:0.01:5;
prior = normpdf(x, -1);
obs = normpdf(x, 1.0);
product = prior .* obs;
plot(x, obs, 'r')
plot(x, prior, 'g')
% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
% Talk about how one might update the ensemble

% Option 1, just form a random sample of the blue curve
% Note that the result is pretty Gaussian (this is NOT EnKF)
% Posterior has mean 1 and variance 1/2 (sd = 1/sqrt(2))
post_x = normrnd(0, 1/sqrt(2), 4, 1)
% Adjust the mean to correct value (0)
post_x = post_x - mean(post_x);
h = plot(post_x, prior_y + 0.05, 'b*');

