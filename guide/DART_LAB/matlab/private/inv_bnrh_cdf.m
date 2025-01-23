function [x] = inv_bnrh_cdf(quantiles, ens_size, sort_ens, ...
   bounded_below, bounded_above, lower_bound, upper_bound, ...
   tail_amp_left, tail_mean_left, tail_sd_left, do_uniform_tail_left, ...
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right)
                 
%% inv_bnrh_cdf Computes the inverse cdf for a bnrh distribution
% This is modified directly from the Fortran version and is not currently
% using matlab efficiently and clearly.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Quantile increment between ensemble members for bnrh
del_q = 1.0 / (ens_size + 1.0);

for i = 1:ens_size
   q(i) = i * del_q;
end

% Loop through each ensemble member to find posterior state
for i = 1:ens_size
   curr_q = quantiles(i);
   % Which region is this quantile in?
   % BNRH quantiles are uniform; finding region for this quantile is trivial
   region = floor(curr_q * (ens_size + 1.0));
   % Careful about numerical issues moving outside of region [0 ens_size]
   if(region < 0) region = 0; end
   if(region > ens_size) region = ens_size; end

   if(region == 0)
      % Lower tail
      if(bounded_below & do_uniform_tail_left)
         % Lower tail uniform
         upper_state = sort_ens(1);
         x(i) = lower_bound + (curr_q / q(1)) * (upper_state - lower_bound);
      else
         % Find the mass at the lower bound (which could be unbounded)
         if(bounded_below)
            lower_mass = tail_amp_left * ...
               normcdf(lower_bound, tail_mean_left, tail_sd_left);
         else
            lower_mass = 0.0;
         end
         % Find the mass at the upper bound (ensemble member 1)
         upper_mass = tail_amp_left * ...
            normcdf(sort_ens(1), tail_mean_left, tail_sd_left);
         % What fraction of this mass difference should we go?
         fract = curr_q / q(1);
         target_mass = lower_mass + fract * (upper_mass - lower_mass);
         x(i) = weighted_norm_inv(tail_amp_left, tail_mean_left, ...
            tail_sd_left, target_mass);
      end

   elseif(region == ens_size)
      % Upper tail
      if(bounded_above & do_uniform_tail_right)
         % Upper tail is uniform
         lower_state = sort_ens(ens_size);
         upper_state = upper_bound;
         x(i) = lower_state + (curr_q - q(ens_size)) * ...
            (upper_state - lower_state) / (1.0 - q(ens_size));
      else
         % Upper tail is (bounded) normal
         % Find the mass at the upper bound (which could be unbounded)
         if(bounded_above) then
            upper_mass = tail_amp_right * ...
               normcdf(upper_bound, tail_mean_right, tail_sd_right);
         else
            upper_mass = 1.0;
         end
         % Find the mass at the lower edge of the region (ensemble member n)
         lower_mass = tail_amp_right * ...
            normcdf(sort_ens(ens_size), tail_mean_right, tail_sd_right);
         % What fraction of the last interval do we need to move
         fract = (curr_q - q(ens_size)) / (1.0 - q(ens_size));
         target_mass = lower_mass + fract * (upper_mass - lower_mass);
         x(i) = weighted_norm_inv(tail_amp_right, tail_mean_right, ...
            tail_sd_right, target_mass);
      end

   else
      % Interior region; get the quantiles of the region boundary
      lower_q = q(region);
      upper_q = q(region + 1);
      x(i) = sort_ens(region) + ((curr_q - lower_q) / (upper_q - lower_q)) * ...
         (sort_ens(region + 1) - sort_ens(region));
   end
  
   % Imprecision can lead to x being slightly out of bounds, fix it to bounds
   % Not needed for matlab test scripting
   %call check_bounds(x(i), curr_q, bounded_below, lower_bound, &
                              %bounded_above, upper_bound, 'inf_bnrh_cdf')
end
