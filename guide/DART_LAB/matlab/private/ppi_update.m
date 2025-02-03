function [post_state, prior_obs_ppi, post_obs_ppi, prior_state_ppi, post_state_ppi] = ...
    ppi_update(prior_obs, prior_state, post_obs, state_dist_type, obs_dist_type)

%% ppi_update Computes increments for unobserved variable with gamma transform
% Need to discuss the available options eventually
% For now this implements the default options

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

ens_size = size(prior_state, 2);

% Convert the prior_obs and posterior_obs to PPI space
% First get quantiles
switch(obs_dist_type)
   case 'Normal'
      % This is just a scaling and removing the mean but done explicitly here
      prior_obs_mean = mean(prior_obs);
      prior_obs_sd = std(prior_obs);
% Note that both are transformed with the prior statistics as per QCEFF algorithm
      prior_obs_q = normcdf(prior_obs, prior_obs_mean, prior_obs_sd);
      post_obs_q = normcdf(post_obs, prior_obs_mean, prior_obs_sd);
   
   case 'RHF'
      % Get the quantiles for this ensemble for an unbounded BNRH distribution
      % Everything except the 
      [sort_x, prior_obs_q, tail_amp_left, tail_mean_left, tail_sd_left, ...
       tail_amp_right, tail_mean_right, tail_sd_right, ...
       do_uniform_tail_left, do_uniform_tail_right] = bnrh_cdf(prior_obs, ...
       ens_size, false, false, -99, -99);

      % Need to use the prior_obs distribution for the quantile transform
      for i = 1:ens_size
         post_obs_q(i) = bnrh_cdf_initialized(post_obs(i), ens_size, sort_x, ...
         false, false, -99, -99, ...
         tail_amp_left, tail_mean_left, tail_sd_left, do_uniform_tail_left, ...
         tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, prior_obs_q);
      end

end

% Do the probit transform
prior_obs_ppi = norminv(prior_obs_q, 0, 1);
post_obs_ppi = norminv(post_obs_q, 0, 1);

% Get the obs_increments in the transformed space
obs_increments = post_obs_ppi - prior_obs_ppi;


% Get the quantiles for the appropriate distribution
switch(state_dist_type)
   case 'Normal'
      % Get the mean and std
      pmean = mean(prior_state);
      psd = std(prior_state);
      % Get quantiles 
      prior_state_q = normcdf(prior_state, pmean, psd);

   case 'Gamma'
      % Get the shape and scale of the prior state ensemble
      pgamma_params = gamfit(prior_state);
      pgamma_shape = pgamma_params(1); pgamma_scale = pgamma_params(2);

      % Compute the quantiles of the prior state members
      prior_state_q = gamcdf(prior_state, pgamma_shape, pgamma_scale);

   case 'RHF'
      % Get the quantiles for this ensemble for an unbounded BNRH distribution
      [sort_x, prior_state_q, tail_amp_left, tail_mean_left, tail_sd_left, ...
       tail_amp_right, tail_mean_right, tail_sd_right, ...
       do_uniform_tail_left, do_uniform_tail_right] = bnrh_cdf(prior_state, ...
       ens_size, false, false, -99, -99);

   case 'BNRH'
      % Get the quantiles for this ensemble for an unbounded BNRH distribution
      [sort_x, prior_state_q, tail_amp_left, tail_mean_left, tail_sd_left, ...
       tail_amp_right, tail_mean_right, tail_sd_right, ...
       do_uniform_tail_left, do_uniform_tail_right] = bnrh_cdf(prior_state, ...
       ens_size, true, false, 0, -99);

end

% Do a probit transform on the quantiles
prior_state_ppi = norminv(prior_state_q, 0, 1);

% Update the unobserved variable
covar = cov(prior_obs_ppi, prior_state_ppi);
state_inc = obs_increments * covar(1, 2) / covar(1, 1);
post_state_ppi = prior_state_ppi + state_inc;

% Transform back to gamma cdf space
post_state_q = normcdf(post_state_ppi, 0, 1);

% Transform back to regular space
switch(state_dist_type)
   case 'Normal'
      post_state = norminv(post_state_q, pmean, psd);

   case 'Gamma'
      post_state = gaminv(post_state_q, pgamma_shape, pgamma_scale);

   case 'RHF'
      post_state = inv_bnrh_cdf(post_state_q, ens_size, sort_x, ...
         false, false, -99, -99, ...
         tail_amp_left, tail_mean_left, tail_sd_left, do_uniform_tail_left, ...
         tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right);

   case 'BNRH'
      post_state = inv_bnrh_cdf(post_state_q, ens_size, sort_x, ...
         true, false, 0, -99, ...
         tail_amp_left, tail_mean_left, tail_sd_left, do_uniform_tail_left, ...
         tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right);
end


%-----------------------------------------------

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
