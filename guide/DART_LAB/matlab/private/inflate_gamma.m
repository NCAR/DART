function [inf_ens] = inflate_gamma(ens, ens_size, var_inf)

%% inflate_gamma Inflates an ensemble in a gamma/probit transformed space

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% The sd inflation
inf = sqrt(var_inf);

% Get the shape and scale of the prior
pgamma_params = gamfit(ens);
pgamma_shape = pgamma_params(1); pgamma_scale = pgamma_params(2);

% Compute the quantiles of the prior ensemble members
prior_q = gamcdf(ens, pgamma_shape, pgamma_scale);

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
inf_ens = gaminv(inf_prior_q, pgamma_shape, pgamma_scale);


