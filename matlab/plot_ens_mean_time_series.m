%% DART:plot_ens_mean_time_series  time series of ensemble mean and truth
%                                                                               
% plot_ens_mean_time_series    interactively queries for the needed information.
%              Since different models potentially need different pieces of 
%              information ... the model types are determined and additional 
%              user input may be queried.
%
% Ultimately, plot_ens_mean_time_series will be replaced by a GUI.
%
%
% All the heavy lifting is done by PlotEnsMeanTimeSeries.

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist('diagn_file','var') ~=1)
   disp(' ')
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

if (exist('truth_file','var') ~= 1)
   disp(' ')
   disp('OPTIONAL: if you have the true state and want it superimposed, provide')
   disp('        : the name of the input file. If not, enter a dummy filename.')
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

vars  = CheckModel(diagn_file);   % also gets default values for this model.

if (exist(truth_file,'file')==2)
   pinfo = CheckModelCompatibility(truth_file, diagn_file);
else
   pinfo.truth_file = [];
end
truth_file  = pinfo.truth_file;


switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection'} 

      varid = SetVariableID(vars);      % queries for variable IDs if needed.
      pinfo = setfield(pinfo, 'truth_file', truth_file);
      pinfo = setfield(pinfo, 'diagn_file', diagn_file);
      pinfo = setfield(pinfo, 'var'       , varid.var);
      pinfo = setfield(pinfo, 'var_inds'  , varid.var_inds);

      fprintf('Comparing %s and \n          %s\n', pinfo.truth_file, pinfo.diagn_file)
      disp(['Using State Variable IDs ', num2str(pinfo.var_inds)])
      clear varid

   case 'fms_bgrid'

      pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');
      pinfo.truth_file = truth_file;   % since it has been verified to be compatible.
      pinfo.diagn_file = diagn_file;   % since it has been verified to be compatible.

   case 'cam'

      vars.truth_file     = truth_file;
      vars.diagn_file     = diagn_file;
      vars.prior_file     = [];
      vars.posterior_file = [];
      pinfo               = GetCamInfo(vars, 'PlotEnsMeanTimeSeries');
   %  pinfo.copyindices   = SetCopyID2(vars.diagn_file);
   %  pinfo.copies        = length(pinfo.copyindices);
      pinfo.truth_file    = truth_file;
      pinfo.diagn_file    = diagn_file;

   case 'wrf'

      pinfo = GetWRFInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;

   case 'pe2lyr'

      pinfo = GetPe2lyrInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;

   case 'mitgcm_ocean'

      pinfo = GetMITgcm_oceanInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');
      pinfo.truth_file = truth_file;
      pinfo.diagn_file = diagn_file;

   otherwise

      error('model %s not implemented yet', vars.model)

end

pinfo

PlotEnsMeanTimeSeries( pinfo )
clear vars
