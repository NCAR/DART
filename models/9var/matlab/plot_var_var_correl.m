% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time  sequence

ens_file = input('Input name of prior or posterior diagnostics file; <cr> for Prior_Diag.nc');
if sum(size(ens_file)) == 0
   ens_file = 'Prior_Diag.nc';
end

base_index = input('Input index for base variable');
base_time = input('Input time index for base point');
var_index = input('Input index for variable to be correlated against');

% Get ensemble mean and spread
base = get_ens_series(base_index, ens_file);
var = get_ens_series(var_index, ens_file);

correl = ens_correl(base, base_time, var);

plot(correl);


