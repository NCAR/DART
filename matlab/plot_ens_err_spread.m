% DART : Plots summary plots of the ensemble error and ensemble spread.
%                        Interactively queries for the needed information.
%                        Since different models potentially need different 
%                        pieces of information ... the model types are 
%                        determined and additional user input may be queried.
%
% Ultimately, plot_ens_err_spread will be replaced by a GUI.
% All the heavy lifting is done by PlotEnsErrSpread.
%
% Example 1 (for low-order models)
%
% truth_file = 'True_State.nc';
% diagn_file = 'Prior_Diag.nc';
% plot_ens_err_spread

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

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_96'}

      pinfo = struct('truth_file'    , truth_file, ...
                     'diagn_file'    , diagn_file, ...
                     'state_var_inds', varid);

      disp(sprintf('Comparing %s and \n          %s', pinfo.truth_file, pinfo.diagn_file))
      disp(['Using State Variable IDs ', num2str(pinfo.state_var_inds)]) 

      PlotEnsErrSpread( pinfo )

   case 'fms_bgrid'

      pinfo = GetBgridInfo(truth_file, 'PlotEnsErrSpread');
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;
      PlotEnsErrSpread( pinfo )

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

clear vars varid
