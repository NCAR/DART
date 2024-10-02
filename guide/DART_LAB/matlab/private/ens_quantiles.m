function [q] = ...
    ens_quantiles(sorted_ens, ens_size, bounded_below, bounded_above, lower_bound, upper_bound)
             

%% ensemble_quantiles Computes the quantiles of a prior ensemble for  a bnrh distribution
% This is modified directly from the Fortran version and is not currently
% using matlab efficiently and clearly.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% This is directly converted from fortran code
% Given sorted ensemble which may have members identical to the bounds or may contain
% duplicates, compute the quantiles for each member in an bounded normal rh distribution

% Get number of ensemble members that are duplicates of the lower bound
lower_dups = 0;
if(bounded_below)
   for i = 1:ens_size
      if(sorted_ens(i) == lower_bound)
         lower_dups = lower_dups + 1;
      else
         break;
      end
   end
end

% Get number of ensemble members that are duplicates of the upper bound
upper_dups = 0;
if(bounded_above)
   for i = ens_size:-1:1
      if(sorted_ens(i) == upper_bound) then
         upper_dups = upper_dups + 1
      else
         break;
      end
   end
end

% If there are duplicate ensemble members away from the boundaries need to revise quantiles
% Make sure not to count duplicates already handled at the boundaries
% Outer loop determines if a series of duplicates starts at sorted index i
d_start = lower_dups + 1;
d_end   = ens_size - upper_dups;

% Get start, length, and end of each series of duplicates away from the bounds
series_num = 1;
series_start(series_num) = d_start;
series_length(series_num) = 1;
for i = d_start + 1:d_end
   if(sorted_ens(i) == sorted_ens(i - 1))
      series_length(series_num) = series_length(series_num) + 1;
   else
      series_end(series_num) = i-1;
      series_num = series_num + 1;
      series_start(series_num) = i;
      series_length(series_num) = 1;
   end
end

% Off the end, finish up the last series
series_end(series_num) = d_end;

% Now get the value of the quantile for the exact ensemble members
% Start with the lower bound duplicates
for i = 1:lower_dups
   q(i) = lower_dups / (2.0 * (ens_size + 1.0));
end

% Top bound duplicates next
for i = ens_size - upper_dups + 1:ens_size
   q(i) = 1.0 - upper_dups / (2.0 * (ens_size + 1.0));
end

% Do the interior series
for i = 1:series_num
   for j = series_start(i):series_end(i)
      q(j) = series_start(i) / (ens_size + 1.0) + (series_length(i) - 1.0) / (2.0 * (ens_size + 1.0));
   end
end
