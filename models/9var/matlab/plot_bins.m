% Plots ensemble rank histograms for 9 variable Lorenz model

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


% Use three different figures with three subplots each
for i = 1:3
   for j = 1:3
      var = (i - 1) + j;
      ens = get_ens_series(var, ens_file);
      truth = get_var_series(var, truth_index, truth_file);
      bins = rank_hist(ens, truth);
      figure(i);
      hold on;
      subplot(3, 1, j);
      bar(bins);
   end
end


