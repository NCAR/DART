function [quantile] = bnrh_cdf_initialized(x, ens_size, sort_ens, ...
   bounded_below, bounded_above, lower_bound, upper_bound, ...
   tail_amp_left, tail_mean_left, tail_sd_left, do_uniform_tail_left, ...
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, q)
             

%% bnrh_cdf_initialed Computes the quantiles for a bnrh distribution for a different 
% ensemble. 
% This is modified directly from the Fortran version and is not currently
% using matlab efficiently and clearly.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Quantile increment between ensemble members for bnrh
del_q = 1.0 / (ens_size + 1.0);

if(x < sort_ens(1)) 
   % In the left tail
   % Do an error check to make sure ensemble member isn't outside bounds, may be redundant
   if(bounded_below & x < lower_bound)
      stop
      %write(errstring, *) 'Ensemble member less than lower bound', x, lower_bound
      %call error_handler(E_ERR, 'bnrh_cdf_initialized', errstring, source)
      % This error can occur due to roundoff in increment generation from BNRHF
      % See discussion in function fix_bounds.
   end

   if(do_uniform_tail_left)
      % Uniform approximation for left tail; Note that denominator cannot be 0 but could be small
      quantile = (x - lower_bound) / (sort_ens(1) - lower_bound) * del_q;
   else
      % It's a normal tail
      if(bounded_below)
         quantile = tail_amp_left * (normcdf(x,           tail_mean_left, tail_sd_left) - ... 
                                     normcdf(lower_bound, tail_mean_left, tail_sd_left));
      else        % Unbounded, tail normal goes all the way down to quantile 0, amplitude is 1
         quantile = (normal_cdf(x,          tail_mean_left, tail_sd_left) / ... 
                    normal_cdf(sort_ens(1), tail_mean_left, tail_sd_left)) * del_q;
      end
      % Make sure it doesn't sneak past the quantile of the smallest ensemble member due to round-off
      quantile = min(quantile, q(1));
   end
elseif(x == sort_ens(1))
   % This takes care of cases where there are multiple bnrh values at the bdry or at first ensemble
   quantile = q(1);
elseif(x > sort_ens(ens_size))
   % In the right tail
   % Do an error check to make sure ensemble member isn't outside bounds, may be redundant
   if(bounded_above & x > upper_bound)
      stop
      %write(errstring, *) 'Ensemble member greater than upper bound first check(see code)', x, upper_bound
      %call error_handler(E_ERR, 'bnrh_cdf_initialized', errstring, source)
      % This error can occur due to roundoff in increment generation from bounded BNRHF
      % See discussion in function fix_bounds
   end

   if(do_uniform_tail_right)
      % Uniform approximation for right tail
      % The division here could be a concern. However, if sort_ens(ens_size) == upper_bound, then
      % x cannot be > sort_ens(ens_size).
      quantile = ens_size *del_q + ... 
         (x - sort_ens(ens_size)) / (upper_bound - sort_ens(ens_size)) * del_q;
   else
      % It's a normal tail
      q_at_largest_ens = normcdf(sort_ens(ens_size), tail_mean_right, tail_sd_right);
      % Want to avoid quantiles exceeding 1 due to numerical issues. Do fraction of the normal part
      if(bounded_above)
         upper_q = tail_amp_right * normcdf(upper_bound, tail_mean_right, tail_sd_right);
         fract = (tail_amp_right * normcdf(x, tail_mean_right, tail_sd_right) - ... 
                  tail_amp_right * q_at_largest_ens) / (upper_q - tail_amp_right * q_at_largest_ens);
      else
         % Normal goes all the way to infinity, amplitude is 1, q at infinity is 1
         fract = (normcdf(x, tail_mean_right, tail_sd_right) - q_at_largest_ens) / ...
            (1.0 -  q_at_largest_ens);
      end

      quantile = ens_size * del_q + fract * del_q;
      quantile = min(quantile, 1.0);
   end

else
   % In an interior bin
   for j = 1:ens_size - 1
      if(x < sort_ens(j+1))
         % The division here could be a concern. 
         % However, sort_ens(j)< x < sort_ens(j+1) so the two cannot be equal
         quantile = j * del_q + ...
            ((x - sort_ens(j)) / (sort_ens(j+1) - sort_ens(j))) * del_q;
         break
      elseif(x == sort_ens(j+1))
         quantile = q(j+1);
         break
      end
   end
end

