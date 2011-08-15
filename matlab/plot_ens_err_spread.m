%% DART: plot_ens_err_spread - summary plots of the ensemble error and ensemble spread.
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

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist('truth_file','var') ~= 1)
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('diagn_file','var') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

pinfo = CheckModelCompatibility(truth_file, diagn_file);
vars  = CheckModel(truth_file);   % also gets default values for this model.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection'}

      varid = SetVariableID(vars);
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;
      pinfo.var        = varid.var;
      pinfo.var_inds   = varid.var_inds;
      %pinfo = struct('truth_file', truth_file, ...
      %               'diagn_file', diagn_file, ...
      %               'var'       , varid.var, ...
      %               'var_inds'  , varid.var_inds);

      fprintf('Comparing %s and \n          %s\n', pinfo.truth_file, pinfo.diagn_file)
      fprintf('Using Variable %s IDs %s\n', pinfo.var,num2str(pinfo.var_inds))
      clear varid

   case 'fms_bgrid'

      pinfo = GetBgridInfo(pinfo, truth_file, 'PlotEnsErrSpread');
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;

   case 'cam'

      pinfo = CombineStructs(pinfo,vars);
      pinfo = GetCamInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case 'pe2lyr'

      pinfo = GetPe2lyrInfo(pinfo, truth_file, 'PlotEnsErrSpread');
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;

   case 'mitgcm_ocean'

      pinfo = GetMITgcm_oceanInfo(pinfo, truth_file, 'PlotEnsErrSpread');
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;

   otherwise

      error('model %s not implemented yet', vars.model)

end

pinfo

PlotEnsErrSpread( pinfo )

clear vars
