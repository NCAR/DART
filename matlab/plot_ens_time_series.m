% DART : Plots time series of ensemble members, mean and truth
%
% Example 2
% truth_file = 'True_State.nc';
% diagn_file = 'Posterior_Diag.nc';
% plot_ens_time_series

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist('truth_file') ~= 1)
   disp('If the True_State.nc exists, it will be plotted. If not, don''t worry.')
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

vars  = CheckModel(diagn_file);   % also gets default values for this model.

if (exist(truth_file)==2)
   pinfo = CheckModelCompatibility(truth_file, diagn_file);
else
   pinfo.truth_file = [];
end
truth_file = pinfo.truth_file;

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'forced_lorenz_96','lorenz_04','ikeda','simple_advection'}

      varid = SetVariableID(vars);      % queries for variable IDs if needed.
      pinfo = setfield(pinfo,'truth_file', truth_file);
      pinfo = setfield(pinfo,'diagn_file', diagn_file);
      pinfo = setfield(pinfo,'var'       , varid.var);
      pinfo = setfield(pinfo,'var_inds'  , varid.var_inds);
      %pinfo = struct('truth_file', truth_file, ...
      %               'diagn_file', diagn_file, ...
      %               'var'       , varid.var, ...
      %               'var_inds'  , varid.var_inds);

    % disp(sprintf('Comparing %s and \n          %s', pinfo.truth_file, pinfo.diagn_file))
    % disp(sprintf('Using Variable %s IDs %s', pinfo.var,num2str(pinfo.var_inds)))

   case 'fms_bgrid'

      varid = SetVariableID(vars);      % queries for variable IDs if needed.
      pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotEnsTimeSeries');
      pinfo.truth_file = truth_file;   % known to be compatible.
      pinfo.diagn_file = diagn_file;   % known to be compatible.

   case 'cam'

      vars.truth_file     = truth_file; 
      vars.diagn_file     = diagn_file; 
      vars.prior_file     = []; 
      vars.posterior_file = []; 
      pinfo               = GetCamInfo(vars, 'PlotEnsTimeSeries');
      pinfo.truth_file    = truth_file;
      pinfo.diagn_file    = diagn_file;

   otherwise
      error(sprintf('model %s not implemented yet', vars.model))

end

pinfo % echo for posterity.

PlotEnsTimeSeries( pinfo )
clear vars varid

