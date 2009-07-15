function corr = ens_correl(base_var, base_time, state_var)
% ens_correl  Computes correlation of a variable at a time to a time series of
% another variable (could be the same one)

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Extract sample of base at base time

base_ens = base_var(base_time, :);

% size(base_var)
% size(base_time)
% size(state_var)
% size(base_ens)

% Loop through time to correlate with the other ensemble series
num_times = size(state_var, 1);
for i = 1:num_times
   x = corrcoef(base_ens, state_var(i, :));
   corr(i) = x(1, 2);
end 
