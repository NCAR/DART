function err = rms_err(pred, verif)
%RMS_ERR: Computes rms error for time series of set of state variables

% Pred and verif are time_series_length x number of variables
num_times = size(pred, 1);
num_vars  = size(pred, 2);

% fprintf('Number of variables in rms err is %d \n', num_vars)

% TJH -- replacing the slow 'loop' 
%err(num_times) = 0.0;
%for i = 1:num_times
%   err(i) = sqrt(sum(abs(pred(i, :) - verif(i, :))));
%end
% with the vectorized version 
err = sqrt( mean( (pred - verif).^2, 2) );
