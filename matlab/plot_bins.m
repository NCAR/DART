% Plots ensemble rank histograms 
%
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Prior_Diag.nc';
% plot_bins

% error_checking ... 
% exist('bob') == 1   means the variable exists. 
%                     the value of the variable is checked later.

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

CheckModelCompatibility(truth_file,diagn_file)
vars  = CheckModel(truth_file);   % also gets default values for this model.
varid = SetVariableID(vars);      % queries for variable IDs if needed.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_96'}

      pinfo = struct('state_var_inds',varid);
      pinfo.truth_file = truth_file;   % since it has been verified to be compatible.
      pinfo.diagn_file = diagn_file;   % since it has been verified to be compatible.

      disp(sprintf('Comparing %s and \n          %s', pinfo.truth_file, pinfo.diagn_file))
      disp(['Using State Variable IDs ', num2str(pinfo.state_var_inds)])

   case 'fms_bgrid'

      pinfo = GetBgridInfo(diagn_file, 'PlotBins');
      pinfo.truth_file = truth_file;   % since it has been verified to be compatible.
      pinfo.diagn_file = diagn_file;   % since it has been verified to be compatible.

      pinfo                            % just echo stuff for posterity.

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

PlotBins(pinfo);
clear vars varid
