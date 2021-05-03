function enh_density = enh_compute_new_density(dist_2, sigma_p_2, sigma_o_2, alpha, beta, gamma_corr, lambda, ens_size)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%

% Compute probability of this lambda being correct
exp_prior = - beta / lambda;

% Compute probability that observation would have been observed given this lambda
fac1 = (1 + gamma_corr * (sqrt(lambda) - 1.0))^2;
fac2 = -1 / ens_size;
if fac1 < abs(fac2), fac2 = 0; end

theta    = sqrt( (fac1+fac2) * sigma_p_2 + sigma_o_2 );
exp_like = - 0.5 * dist_2 / theta^2;

% Compute the updated probability density for lambda
enh_density = beta^alpha / gamma(alpha) * lambda^(- alpha - 1) / ... 
              (sqrt(2.0 * pi) * theta) * exp(exp_like + exp_prior);
          
