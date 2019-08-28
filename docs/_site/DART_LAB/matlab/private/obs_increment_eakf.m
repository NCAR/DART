function [obs_increments, err] =  obs_increment_eakf(ensemble, observation, obs_error_var)
%% obs_increment_eakf Computes increments for an ensemble adjustment filter

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Set error return to default successful
err = 0;

% Compute prior ensemble mean and variance
prior_mean = mean(ensemble);
prior_var = var(ensemble);

% If both prior and observation error variance are return error
if(prior_var <= 0 && obs_error_var <= 0)
   err = 1;
   return;
end

% Compute the posterior mean and variance
% If prior variance is 0, posterior mean is prior_mean and variance is 0
if(prior_var == 0)
   post_mean = prior_mean;
   post_var = 0;
elseif(obs_error_var == 0)
% If obs_error_var is 0, posterior mean is observation and variance is 0
   post_mean = observation;
   post_var = 0;
else
% Use product of gaussians
   % Compute the posterior variance
   post_var = 1 / (1 / prior_var + 1 / obs_error_var);

   % Compute posterior mean
   post_mean = post_var * (prior_mean / prior_var + observation / obs_error_var);
end

% Shift the prior ensemble to have the posterior mean
updated_ensemble = ensemble - prior_mean + post_mean;

% Contract the ensemble to have the posterior_variance
var_ratio = post_var / prior_var;
updated_ensemble = sqrt(var_ratio) * (updated_ensemble - post_mean) + post_mean;

% Compute the increments
obs_increments = updated_ensemble - ensemble;


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
