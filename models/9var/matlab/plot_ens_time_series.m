% Plots summary plots of error and spread for Lorenz 9 variable model

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

%err = rms_err(ens, truth);

spread = get_state_copy(ens_spread_index, ens_file);

% Use three different figures with three subplots each
for i = 1:3
   for j = 1:3
      var = (i - 1)*3 + j;
      err = rms_err(ens(:, var), truth(:, var));
      figure(i);
      hold on;
      subplot(3, 1, j);
      plot(err, 'b');
      hold on;
      plot(spread(:, var), 'r');
   end
end


