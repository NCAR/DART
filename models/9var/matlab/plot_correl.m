% Plots space-time series of correlation between a given variable at a given
% time and all other variable at all times in an ensemble time  sequence

ens_file = input('Input name of prior or posterior diagnostics file; <cr> for Prior_Diag.nc');
if sum(size(ens_file)) == 0
   ens_file = 'Prior_Diag.nc';
end

base_index = input('Input index for base variable');
base_time = input('Input time index for base point');

% Get ensemble mean and spread
base = get_ens_series(base_index, ens_file);


% Figure out how many variables there are all together
num_vars = 9

% Need to loop through each other variable in the ensemble
for i = 1:9
   var = get_ens_series(i, ens_file);
   correl(i, :) = ens_correl(base, base_time, var);
end


contour(correl);


