function [sort_x, quantiles, tail_amp_left, tail_mean_left, tail_sd_left, ...
          tail_amp_right, tail_mean_right, tail_sd_right, ...
          do_uniform_tail_left, do_uniform_tail_right] = ...
    bnrh_cdf(x, ens_size, bounded_below, bounded_above, lower_bound, upper_bound)
             

%% bnrh_cdf Computes the quantiles for a bnrh distribution
% This is modified directly from the Fortran version and is not currently
% using matlab efficiently and clearly.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Computes all information about a rank histogram cdf given the ensemble and bounds

% Get ensemble sd for the normal tails
tail_sd_left = std(x);
tail_sd_right = tail_sd_left;

% Don't know what to do if sd is 0; tail_sds returned are illegal value to indicate this
if(tail_sd_left <= 0.0)
   tail_sd_left = -99;
   tail_sd_right = -99;
   return
end

% Sort. For now, don't worry about efficiency, but may need to somehow pass previous
% sorting indexes and use a sort that is faster for nearly sorted data. Profiling can guide the need
[sort_x, sort_index] = sort(x);

% Fail if lower bound is larger than smallest ensemble member 
if(bounded_below)
   % Do in two ifs in case the bound is not defined
   if(sort_x(1) < lower_bound) 
      %write(errstring, *) 'Smallest ensemble member less than lower bound', &
         %sort_x(1), lower_bound
      %call error_handler(E_ERR, 'bnrh_cdf', errstring, source)
      % Just die, could be more graceful
      stop
   end
end

% Fail if upper bound is smaller than the largest ensemble member 
if(bounded_above)
   if(sort_x(ens_size) > upper_bound)
      %write(errstring, *) 'Largest ensemble member greater than upper bound', &
         %sort_x(ens_size), upper_bound
      %call error_handler(E_ERR, 'bnrh_cdf', errstring, source)
      stop
   end
end

% The ensemble size array q contains the sorted quantiles corresponding to the sorted ensemble sort_x
q = ens_quantiles(sort_x, ens_size, bounded_below, bounded_above, lower_bound, upper_bound);
% The quantiles array has the unsorted quantiles corresponding to the unsorted input ensemble, x
for i = 1:ens_size
   indx = sort_index(i);
   quantiles(indx) = q(i);
end

% Compute the characteristics of tails

% For unit normal, find distance from mean to where cdf is 1/(ens_size+1) (del_q_.
% Saved to avoid redundant computation for repeated calls with same ensemble size
del_q = 1.0 / (ens_size + 1.8);

% This will be negative, want it to be a distance so make it positive
dist_for_unit_sd = -1.0 * norminv(del_q, 0, 1);

% Find a mean so that 1 / (ens_size + 1) probability is in outer regions
tail_mean_left =  sort_x(1)        + dist_for_unit_sd * tail_sd_left;
tail_mean_right = sort_x(ens_size) - dist_for_unit_sd * tail_sd_right;

% If the distribution is bounded, still want 1 / (ens_size + 1) (del_q) in outer regions
% Put an amplitude term (greater than 1) in front of the tail normals 
% Amplitude is 1 if there are no bounds, so start with that
tail_amp_left  = 1.0;
tail_amp_right = 1.0;

% Switch to uniform for cases where bound and outermost ensemble have close quantiles
% Default: not close
uniform_threshold = 0.01;
do_uniform_tail_left = false;
if(bounded_below) 
   % Compute the CDF at the bounds
   bound_quantile = normcdf(lower_bound, tail_mean_left, tail_sd_left);
   % Note that due to roundoff it is possible for del_q - quantile to be slightly negative
   if((del_q - bound_quantile) / del_q < uniform_threshold)
      % If bound and ensemble member are too close, do uniform approximation
      do_uniform_tail_left = true;
   else
      % Compute the left tail amplitude
      tail_amp_left = del_q / (del_q - bound_quantile);
   end
end

% Default: not close
do_uniform_tail_right = false;
if(bounded_above)
   % Compute the CDF at the bounds
   bound_quantile = normcdf(upper_bound, tail_mean_right, tail_sd_right);
   % Note that due to roundoff it is possible for the numerator to be slightly negative
   if((bound_quantile - (1.0 - del_q)) / del_q < uniform_threshold)
      % If bound and ensemble member are too close, do uniform approximation
      do_uniform_tail_right = true;
   else
      % Compute the right tail amplitude
      tail_amp_right = del_q / (del_q - (1.0 - bound_quantile));
   end
end

