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

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist('diagn_file','var') ~=1)
   disp(' ')
   disp('Input name of prior or posterior diagnostics file:')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

if (exist('truth_file','var') ~= 1)
   disp(' ')
   disp('OPTIONAL: if you have the true state and want it superimposed, provide')
   disp('        : the name of the input file. If not, enter a dummy filename.')
   disp('        : Input name of True State file:')
   truth_file = input('<cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

pinfo  = CheckModel(diagn_file); % also gets default values for this model.

if (exist(truth_file,'file')==2)

   pinfo  = rmfield(pinfo,{'time_series_length','time','fname'});
   vars   = CheckModelCompatibility(truth_file, diagn_file);
   pinfo  = CombineStructs(pinfo,vars);
   clear vars

else
   truth_file = [];
end

pinfo.diagn_file = diagn_file;
pinfo.truth_file = truth_file;

clear MyInfo mynames myname ifield

%% For each model, do what needs to be done.

switch lower(pinfo.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'forced_lorenz_96','lorenz_04','ikeda','simple_advection'}

      varid = SetVariableID(pinfo);      % queries for variable IDs if needed.
      pinfo.var        = varid.var;
      pinfo.var_inds   = varid.var_inds;

      fprintf('Comparing %s and \n          %s\n', pinfo.truth_file, pinfo.diagn_file)
      disp(['Using State Variable IDs ', num2str(pinfo.var_inds)])
      clear varid

   case {'fms_bgrid'}

      pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   case {'cam'}

      pinfo = GetCamInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   case {'wrf'}

      pinfo = GetWRFInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   case {'pe2lyr'}

      pinfo = GetPe2lyrInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   case {'mitgcm_ocean'}

      pinfo = GetMITgcm_oceanInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   case {'mpas_atm'}

      pinfo = GetMPAS_ATMInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   case {'sqg'}

      pinfo = GetSqgInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   case {'pop'}

      pinfo = GetPOPInfo(pinfo, diagn_file, 'PlotEnsMeanTimeSeries');

   otherwise

      error('model %s not implemented yet', pinfo.model)

end

pinfo

PlotEnsMeanTimeSeries( pinfo )


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

