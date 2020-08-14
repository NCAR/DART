function [new_cov_inflate, new_cov_inflate_sd] = update_inflate(x, sigma_p_2, obs, sigma_o_2, inflate_prior_val, lambda_mean, ...
                                                    lambda_mean_LB, lambda_mean_UB, gamma, lambda_sd, lambda_sd_LB)
% Adaptive scheme to update the inflation both in space and time. 
% Based on the algorithm in the following article:
% Anderson, J. L., 2009: Spatially and temporally varying adaptive covariance 
% inflation for ensemble filters. Tellus A, 61, 72-83. doi: 10.1111/j.1600-0870.2008.00361.x 
%
% More documentation available at:
% https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/assimilation_code/modules/assimilation/adaptive_inflate_mod.html

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$


%% FIRST, update the inflation mean:

% Get the "non-inflated" variance of the sample
% lambda here, is the prior value before the update.
sigma_p_2 = sigma_p_2 / ( 1 + gamma*(sqrt(inflate_prior_val) - 1) )^2;

% Squared-innovation
dist_2 = (x - obs)^2;

% d (distance between obs and ens) is drawn from Gaussian with 0 mean and
% variance: E(d^2) = \lambda^o*\sigma_p^2 + \sigma_o^2
theta_bar_2  = ( 1 + gamma * (sqrt(lambda_mean) - 1) )^2 * sigma_p_2 + sigma_o_2;
theta_bar    = sqrt(theta_bar_2);
u_bar        = 1 / (sqrt(2 * pi) * theta_bar);
like_exp_bar = - 0.5 * dist_2 / theta_bar_2;
v_bar        = exp(like_exp_bar);

gamma_terms  = 1 - gamma + gamma*sqrt(lambda_mean);
dtheta_dinf  = 0.5 * sigma_p_2 * gamma * gamma_terms / (theta_bar * sqrt(lambda_mean));              

% The likelihood: p(d/lambda)
like_bar     = u_bar * v_bar;

% Derivative of the likelihood; evaluated at the current inflation mean
like_prime   = (like_bar * dtheta_dinf / theta_bar) * (dist_2 / theta_bar_2 - 1);
like_ratio   = like_bar / like_prime;

% Solve a quadratic equation
a = 1;
b = like_ratio - 2*lambda_mean;
c = lambda_mean^2 - lambda_sd^2 - like_ratio*lambda_mean ;

o = max( [ abs(a), abs(b), abs(c) ] );
a = a/o;
b = b/o;
c = c/o;
d = b^2 - 4*a*c;

if b < 0
    s1 = 0.5 * ( -b + sqrt(d) ) / a;
else 
    s1 = 0.5 * ( -b - sqrt(d) ) / a;
end
s2 = ( c/a ) / s1;

% Select the updated-mean that is closest to the prior
if abs(s2 - lambda_mean) < abs(s1 - lambda_mean)
   new_cov_inflate = s2;
else
   new_cov_inflate = s1;
end

% Make sure the update is not smaller than the lower bound
if new_cov_inflate < lambda_mean_LB || new_cov_inflate > lambda_mean_UB || isnan(new_cov_inflate)
    new_cov_inflate = lambda_mean_LB; 
    new_cov_inflate_sd = lambda_sd;
    return
end


%% Now, update the inflation variance
if lambda_sd <= lambda_sd_LB 
    new_cov_inflate_sd = lambda_sd;
    return
else
    % First compute the new_max value for normalization purposes
    new_max = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, new_cov_inflate);

    % Find value at a point one OLD sd above new mean value
    new_1_sd = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, new_cov_inflate + lambda_sd);

    ratio = new_1_sd / new_max;

    % Can now compute the standard deviation consistent with this as
    % sigma = sqrt(-x^2 / (2 ln(r))  where r is ratio and x is lambda_sd (distance from mean)
    new_cov_inflate_sd = sqrt( - 0.5 * lambda_sd^2 / log(ratio) );

    % Prevent an increase in the sd of lambda
    if new_cov_inflate_sd > lambda_sd, new_cov_inflate_sd = lambda_sd; end
    if new_cov_inflate_sd < lambda_sd_LB, new_cov_inflate_sd = lambda_sd_LB; end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
