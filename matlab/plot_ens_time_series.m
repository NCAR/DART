% Plots time series of ensemble members, mean and truth
%
% Example 2
% truth_file = 'True_State.nc';
% diagn_file = 'Posterior_Diag.nc';
% plot_ens_time_series

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

CheckModelCompatibility(truth_file, diagn_file)
vars  = CheckModel(truth_file);   % also gets default values for this model.
varid = SetVariableID(vars);      % queries for variable IDs if needed.

disp(sprintf('Comparing %s and \n          %s',truth_file,diagn_file))
disp(['Using State Variable IDs ',num2str(varid)])

PlotEnsTimeSeries(truth_file, diagn_file, varid)
clear vars varid

