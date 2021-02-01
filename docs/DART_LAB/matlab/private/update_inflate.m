function [new_cov_inflate, new_cov_inflate_sd] = update_inflate(x_p, r_var, y_o, sigma_o_2, ss_inflate_base, lambda_mean, lambda_sd, ...
                                                    inf_lower_bound, inf_upper_bound, gamma_corr, sd_lower_bound_in, ens_size, flavor)
% Adaptive scheme to update the inflation both in space and time. 
% Based on the algorithm in the following articles:
%
% [flavor: 2, Gaussian inflation pdf] 
% Anderson, J. L., 2009: Spatially and temporally varying adaptive covariance 
% inflation for ensemble filters. Tellus A, 61, 72-83. doi: 10.1111/j.1600-0870.2008.00361.x 
%
% [flavor: 5, I-Gamma inflation pdf] 
% El Gharamti, M., 2018: Enhanced Adaptive Inflation Algorithm for Ensemble Filters. 
% Monthly Weather Review, 146, 623-640. doi: 10.1175/MWR-D-17-0187.1
%
% More documentation available at:
% assimilation_code/modules/assimilation/adaptive_inflate_mod.html

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

%% FIRST, update the inflation mean:

% Get the "non-inflated" variance of the sample
% lambda here, is the prior value before the update.
sigma_p_2 = r_var / ( 1 + gamma_corr*(sqrt(ss_inflate_base) - 1) )^2;

% Squared-innovation
dist_2 = (x_p - y_o)^2;

% d (distance between obs and ens) is drawn from Gaussian with 0 mean and
% variance: E(d^2) = \lambda^o*\sigma_p^2 + \sigma_o^2
if strcmp(flavor, 'Gaussian')
    theta_bar_2  = ( 1 + gamma_corr * (sqrt(lambda_mean) - 1) )^2 * sigma_p_2 + sigma_o_2;
    
elseif strcmp(flavor, 'I-Gamma')
    fac1 = (1 + gamma_corr * (sqrt(lambda_mean) - 1))^2;
    fac2 = -1 / ens_size;
    
    if fac1 < abs(fac2), fac2 = 0; end
    
    theta_bar_2 = (fac1+fac2) * sigma_p_2 + sigma_o_2;
end

theta_bar    = sqrt(theta_bar_2);
u_bar        = 1 / (sqrt(2 * pi) * theta_bar);
like_exp_bar = - 0.5 * dist_2 / theta_bar_2;
v_bar        = exp(like_exp_bar);            

% The likelihood: p(d/lambda)
like_bar     = u_bar * v_bar;

% likelihood can't be less than or equal to zero!
if like_bar <= 0
    new_cov_inflate     = inf_lower_bound; 
    new_cov_inflate_sd  = lambda_sd;
    return
end

gamma_terms     = 1 - gamma_corr + gamma_corr*sqrt(lambda_mean);
dtheta_dlambda  = 0.5 * sigma_p_2 * gamma_corr * gamma_terms / (theta_bar * sqrt(lambda_mean));  

% Derivative of the likelihood; evaluated at the current inflation mean
like_prime   = (like_bar * dtheta_dlambda / theta_bar) * (dist_2 / theta_bar_2 - 1);

% Make sure we don't divide by zero
if like_prime == 0
    new_cov_inflate     = inf_lower_bound; 
    new_cov_inflate_sd  = lambda_sd;
    return
end

like_ratio   = like_bar / like_prime;

% Solve a quadratic equation
if strcmp(flavor, 'Gaussian')
    a = 1;
    b = like_ratio - 2*lambda_mean;
    c = lambda_mean^2 - lambda_sd^2 - like_ratio*lambda_mean;
    
elseif strcmp(flavor, 'I-Gamma')
    % Switch from Gaussian prior to inverse-gamma
    beta = change_GA_IG(lambda_mean, lambda_sd^2);
    
    a = 1 - lambda_mean / beta;
    b = like_ratio - 2 * lambda_mean;
    c = lambda_mean^2 - like_ratio * lambda_mean;
end

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
if new_cov_inflate < inf_lower_bound || new_cov_inflate > inf_upper_bound || isnan(new_cov_inflate)
    new_cov_inflate     = inf_lower_bound; 
    new_cov_inflate_sd  = lambda_sd;
    return
end


%% SECOND, update the inflation variance
if lambda_sd <= sd_lower_bound_in 
    new_cov_inflate_sd = lambda_sd;
    return
else
    if strcmp(flavor, 'Gaussian')
        % First compute the new_max value for normalization purposes
        new_max = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma_corr, new_cov_inflate);

        % Find value at a point one OLD sd above new mean value
        new_1_sd = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma_corr, new_cov_inflate + lambda_sd);
        
        % If either the numerator or denominator of the following computation 
        % of 'ratio' is going to be zero (or almost so), return the original incoming
        % inflation value.  The computation would have resulted in either Inf or NaN.
        if abs(new_max) <= realmin || abs(new_1_sd) <= realmin
            new_cov_inflate_sd = lambda_sd;
            return
        end

        ratio = new_1_sd / new_max;
        
        % Another error for numerical issues; if ratio is larger than 0.99, bail out
        if ratio > 0.99 
            new_cov_inflate_sd = lambda_sd;
            return
        end

        % Can now compute the standard deviation consistent with this as
        % sigma = sqrt(-x^2 / (2 ln(r))  where r is ratio and x is lambda_sd (distance from mean)
        new_cov_inflate_sd = sqrt( - 0.5 * lambda_sd^2 / log(ratio) );

        % Prevent an increase in the sd of lambda
        if new_cov_inflate_sd > lambda_sd,          new_cov_inflate_sd = lambda_sd;         end
        
    elseif strcmp(flavor, 'I-Gamma')
        % Compute the shape parameter of the prior IG
        % This comes from the assumption that the mode of the IG is the mean/mode of the input Gaussian
        shape_old = beta / lambda_mean - 1;
        if shape_old <= 2
            new_cov_inflate_sd = lambda_sd;
            return
        end 
        
        % Evaluate the exact IG posterior at p1: \lambda_u+\sigma_{\lambda_b} & p2: \lambda_u
        density_1 = enh_compute_new_density(dist_2, sigma_p_2, sigma_o_2, shape_old, ... 
                                            beta, gamma_corr, new_cov_inflate+lambda_sd, ens_size);
        density_2 = enh_compute_new_density(dist_2, sigma_p_2, sigma_o_2, shape_old, ... 
                                            beta, gamma_corr, new_cov_inflate, ens_size);

        % Computational errors check (small numbers + NaNs)
        if abs(density_1) <= realmin || abs(density_2) <= realmin || isnan(density_1) || isnan(density_2)
             new_cov_inflate_sd = lambda_sd;
             return
        end 
        
        % Now, compute omega and the new distribution parameters
        ratio     = density_1 / density_2;
        omega     = log(new_cov_inflate          )/new_cov_inflate + 1/new_cov_inflate - ...
                    log(new_cov_inflate+lambda_sd)/new_cov_inflate - 1/(new_cov_inflate+lambda_sd);
        rate_new  = log(ratio) / omega;
        shape_new = rate_new / new_cov_inflate - 1;
        
        % Finally, get the sd of the IG posterior
        if (shape_new <= 2)
            new_cov_inflate_sd = lambda_sd;
            return
        end
        new_cov_inflate_sd = sqrt(rate_new^2 / ( (shape_new-1)^2 * (shape_new-2) ));
        
        % If the updated variance is more than 5% the prior variance, keep the prior unchanged 
        % for stability reasons. Also, if the updated variance is NaN (not sure why this
        % can happen; never did when develping this code), keep the prior variance unchanged. 
        if ( new_cov_inflate_sd > 1.05*lambda_sd || isnan(new_cov_inflate_sd) )
            new_cov_inflate_sd = lambda_sd;
        end 
    end
    
    if new_cov_inflate_sd < sd_lower_bound_in,  new_cov_inflate_sd = sd_lower_bound_in; end
end

