% Plots ensemble rank histograms for 9 variable Lorenz model
%
% Example 1
% truth_file = '../work/True_State.nc';
% diagn_file = '../work/Prior_Diag.nc';
% plot_bins

% error_checking ... 
% exist('bob') == 1   means the variable exists. 
%                     the value of the variable is checked later.

if (exist('truth_file') ~= 1) 
   truth_file = input('Input name of True State file; <cr> for ../work/True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = '../work/True_State.nc';
   end
end

if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;') 
   diagn_file = input('<cr> for ../work/Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = '../work/Prior_Diag.nc';
   end
end

disp(sprintf('Comparing %s and \n          %s',truth_file,diagn_file))

PlotBins(truth_file,diagn_file);
