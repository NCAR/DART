% Computes correlation of a variable at a time to a time series of
% another variable (could be the same one)

function corr = ens_correl(base_var, base_time, var)

%Extract sample of base at base time
base_ens = base_var(base_time, :);

% Loop through time to correlate with the other ensemble series
num_times = size(var, 1);
for i = 1:num_times
   x = corrcoef(base_ens, var(i, :));
   corr(i) = x(1, 2);
end 

