% Plots summary plots of error and spread 

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

vars = CheckModel(truth_file);   % also gets default values for this file.

disp(sprintf('Comparing %s and \n          %s',truth_file,diagn_file))
disp(['using variables ',num2str(vars.def_state_vars)]) 

% here's where the "if you don't like these values, input something else"
% block should go ...
% vars.def_state_vars = [ 1 2 3 34 35 18];

PlotEnsErrSpread(truth_file,diagn_file, vars.def_state_vars)

clear vars
