function [post_state, prior_obs_ppi, post_obs_ppi, prior_state_ppi, post_state_ppi] = ...
    gamma_ppi_update(prior_obs, prior_state, post_obs)

%% gamma_ppi_update Computes increments for unobserved variable with gamma transform
% Need to discuss the available options eventually
% For now this implements the default options

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Do the transform for the normally distributed observed ensemble (just a normalizaion)
prior_obs_mean = mean(prior_obs);
prior_obs_sd = std(prior_obs);
% Note that both are transformed with the prior statistics as per QCEFF algorithm
prior_obs_ppi = (prior_obs - prior_obs_mean) / prior_obs_sd + prior_obs_mean;
post_obs_ppi = (post_obs - prior_obs_mean) / prior_obs_sd + prior_obs_mean;

% Get the obs_increments in the transformed space
obs_increments = post_obs_ppi - prior_obs_ppi;

% Get the shape and scale of the prior state ensemble
pgamma_params = gamfit(prior_state);
pgamma_shape = pgamma_params(1); pgamma_scale = pgamma_params(2);

% Compute the quantiles of the prior state members
prior_state_q = gamcdf(prior_state, pgamma_shape, pgamma_scale);

% Do a probit transform on the quantiles
prior_state_ppi = norminv(prior_state_q, 0, 1);

% Update the unobserved variable
covar = cov(prior_obs_ppi, prior_state_ppi);
state_inc = obs_increments * covar(1, 2) / covar(1, 1);
post_state_ppi = prior_state_ppi + state_inc;

% Transform back to gamma cdf space
post_state_q = normcdf(post_state_ppi, 0, 1);

% Transform back to regular space
post_state = gaminv(post_state_q, pgamma_shape, pgamma_scale);

%-----------------------------------------------

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
