function err = total_err(pred, verif)
%TOTAL_ERR: Computes Total error for time series of set of state variables

% Pred and verif are time_series_length x number of variables
num_times = size(pred, 1);
num_vars  = size(pred, 2);

err = sqrt( sum( (pred - verif).^2, 2) );
