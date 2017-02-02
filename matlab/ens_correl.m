function corr = ens_correl(base_var, base_time, state_var)
%% ens_correl  Computes correlation of a variable at a time to a time series of
% another variable (could be the same one)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Extract sample of base at base time

base_ens  = base_var(base_time, :);

% preallocate space for result

num_times = size(state_var, 1);
corr      = zeros(num_times,1);
corr(:)   = NaN;

% Loop through time to correlate with the other ensemble series
for i = 1:num_times
   x = corrcoef(base_ens, state_var(i, :));
   corr(i) = x(1, 2);
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
