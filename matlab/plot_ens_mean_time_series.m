% plot_ens_mean_time_series.m
% Plots time series of ensemble members, mean and truth for 9 variable

if (exist('truth_file') ~= 1)
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

disp(sprintf('Comparing %s and \n          %s',truth_file,diagn_file))

PlotEnsMeanTimeSeries(truth_file,diagn_file)
