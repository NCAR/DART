function [post_mean post_sd weight] =...
   product_of_gaussians(prior_mean, prior_sd, obs, obs_err_sd) 
% Computes mean, variance and weight of the product of two unit gaussians given the mean and standard deviation of each.

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Get the prior and observational error variance
prior_var = prior_sd^2;
obs_err_var = obs_err_sd^2;

% Compute the posterior variance
post_var = 1. / (1. / prior_var + 1. / obs_err_var);
post_sd = sqrt(post_var);

% Compute the posterior mean
post_mean = post_var * (prior_mean / prior_var + obs / obs_err_var);

% Compute the associated weight
weight = (1. / (sqrt(2. * pi) * sqrt(prior_var + obs_err_var))) *...
   exp(-0.5 * (obs - prior_mean).^2 ./ (prior_var + obs_err_var));

