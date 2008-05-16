% DART : Plots ensemble rank histograms 
%
% plot_bins    interactively queries for the information needed to create
%              ensemble rank histograms. Since different models potentially
%              need different pieces of information ... the model types are
%              determined and additional user input may be queried.
%
% Ultimately, plot_bins will be replaced by a GUI.
% In the end, the heavy lifting is done by PlotBins.
%
% Example 1 (for low-order models)
%
% truth_file = 'True_State.nc';
% diagn_file = 'Prior_Diag.nc';
% plot_bins

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

pinfo = CheckModelCompatibility(truth_file,diagn_file);
vars  = CheckModel(truth_file);   % also gets default values for this model.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04', 'forced_lorenz_96','ikeda','simple_advection'}

      varid = SetVariableID(vars);      % queries for variable IDs
      pinfo = setfield(pinfo, 'truth_file'    , truth_file);
      pinfo = setfield(pinfo, 'diagn_file'    , diagn_file);
      pinfo = setfield(pinfo, 'var'           , varid.var);
      pinfo = setfield(pinfo, 'state_var_inds', varid.var_inds);
      clear varid

   case 'fms_bgrid'

      pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotBins');

   case 'cam'

      pinfo = GetCamInfo(pinfo, diagn_file, 'PlotBins');

   case 'pe2lyr'

      pinfo = GetPe2lyrInfo(pinfo, diagn_file, 'PlotBins');

   case 'mitgcm_ocean'

      pinfo = GetMITgcm_oceanInfo(pinfo, diagn_file, 'PlotBins');

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

PlotBins(pinfo);
clear vars
