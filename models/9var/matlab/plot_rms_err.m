% Plots summary plots of global error and spread

truth_file = input('Input name of True State file; <cr> for True_State.nc');
if sum(size(truth_file)) == 0
   truth_file = 'True_State.nc';
end

ens_file = input('Input name of prior or posterior diagnostics file; <cr> for Prior_Diag.nc');
if sum(size(ens_file)) == 0
   ens_file = 'Prior_Diag.nc';
end


% Get the state for the truth
truth_index = get_copy_index('true state', truth_file);

% Get ensemble mean and spread
ens_mean_index = get_copy_index('ensemble mean', ens_file);
ens_spread_index = get_copy_index('ensemble spread', ens_file);

% Get the truth and ens_mean and compute error
truth = get_state_copy(truth_index, truth_file);
ens = get_state_copy(ens_mean_index, ens_file);

spread = get_state_copy(ens_spread_index, ens_file);

% Compute RMS
rms = rms_err(truth, ens);
plot(rms, 'b');

% Also need to compute the spread; zero truth for this and compute
% RMS distance from 0 (verify algorithm)
truth(:) = 0.0;
rms_spread = rms_err(truth, spread);
hold on;
plot(rms_spread, 'r');



