function density = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, lambda)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$


% Compute probability of this lambda being correct
exponent_prior = - 0.5 * (lambda - lambda_mean)^2 / lambda_sd^2;

% Compute probability that observation would have been observed given this lambda
theta_2 = (1.0 + gamma * (sqrt(lambda) - 1.0))^2 * sigma_p_2 + sigma_o_2;
theta = sqrt(theta_2);

exponent_likelihood = dist_2 / ( -2.0 * theta_2);

% Compute the updated probability density for lambda
% Have 1 / sqrt(2 PI) twice, so product is 1 / (2 PI)
density = exp(exponent_likelihood + exponent_prior) / (2.0 * pi * lambda_sd * theta);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
