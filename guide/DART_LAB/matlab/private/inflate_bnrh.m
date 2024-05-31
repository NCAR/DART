function [inf_ens] = inflate_bnrh(ens, ens_size, var_inf, ...
   bounded_below, bounded_above, lower_bound, upper_bound)

%% inflate_bnrh Inflates an ensemble in a bnrh/probit transformed space

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% The sd inflation
inf = sqrt(var_inf);

% Compute the quantiles of the prior ensemble members
[sort_x, prior_q, tail_amp_left, tail_mean_left, tail_sd_left, ...
   tail_amp_right, tail_mean_right, tail_sd_right, ...
   do_uniform_tail_left, do_uniform_tail_right] = ...
   bnrh_cdf(ens, ens_size, bounded_below, bounded_above, lower_bound, upper_bound);

% Do a probit transform on the quantiles
probit_ens = norminv(prior_q, 0, 1);

% Inflate in probit space
probit_inf_ens = zeros(1, ens_size);
probit_mean = mean(probit_ens);
for i = 1:ens_size
   probit_inf_ens(i) = (probit_ens(i) - probit_mean) * inf + probit_mean;
end

% Transform back to gamma cdf space
inf_prior_q = normcdf(probit_inf_ens, 0, 1);

% Transform back to regular space
inf_ens = inv_bnrh_cdf(inf_prior_q, ens_size, sort_x, ...
   bounded_below, bounded_above, lower_bound, upper_bound, ...
   tail_amp_left, tail_mean_left, tail_sd_left, do_uniform_tail_left, ...
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right);

