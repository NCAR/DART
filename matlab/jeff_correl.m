function corr = ens_correl(base_var, base_time, state_var)
%% ens_correl  Computes correlation of a variable at a time to a time series of
% another variable (could be the same one)

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

%Extract sample of base at base time
%base_ens = base_var(base_time, :);

% size(base_var)
% size(base_time)
% size(state_var)
% size(base_ens)

% Loop through time to correlate with the other ensemble series
figure(9);
hold on;
num_times = size(state_var, 1);
for i = 1:num_times
   x = corrcoef(base_var(i, :), state_var(i, :));
   corr(i) = x(1, 2);
%   plot(base_var(i, 1), state_var(i, 1), 'x');
%   plot(base_var(i, :), state_var(i, :), 'x');
%   plot(2, corr(i), 'x');
end 

plot(base_var(:, 1), state_var(:, 1));

figure(1);
