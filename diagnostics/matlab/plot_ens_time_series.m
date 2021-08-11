function plot_ens_time_series(input_file)

%% DART:plot_ens_time_series - time series of ensemble and truth (if available)
%
% plot_ens_time_series    interactively queries for the needed information.
%              Since different models potentially need different pieces of
%              information ... the model types are determined and additional
%              user input may be queried.
%              The first filename specifies the file with the ensemble mean.
%              If no filename is provided, you will be prompted for one.
%              If the 'true_state.nc' is available, it is used to plot the truth.
%
% The true state file is specified by the variable 'truth_file'.
%
% A reminder of the sequence:
% truth  run (from    pmo):
%           perfect_input.nc --->     true_state.nc     ---> perfect_output.nc
%           [single timestep]      [multiple timesteps]      [single timestep]
%
% filter run (from filter):
%           filter_input.nc  --->
%               forecast.nc      ---> [prior inflation] --->
%               ^     preassim.nc     ---> [assimilation] --->
%               ^          postassim.nc    ---> [posterior inflation] --->
%               ^<<<<<<<<<<<<<<< analysis.nc   --->
%                                      filter_output.nc
%
% Example 1:  Prompt for filter output filename. Default is defined in 'diagn_file' variable.
%             The true state filename always comes from the 'truth_file' variable.
% plot_ens_time_series
%
% Example 2:  Prompt, but use name specified in 'diagn_file'.
% diagn_file = 'analysis.nc';
% plot_ens_time_series
%
% Example 3:  Prompt with 'diagn_file', different true state.
% diagn_file = 'analysis.nc';
% truth_file = 'More_True_State.nc';
% plot_ens_time_series
%
% Example 4: No prompting, true state filename comes from 'truth_file' variable.
% plot_ens_time_series('output.nc')

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

[truth_file, diagn_file] = set_default_files();

if (nargin == 0)
    mystring = sprintf('<cr> for %s\n',diagn_file);
    disp('Input name of ensemble trajectory file:')
    input_file = input(mystring,'s');
    if isempty(input_file)
        input_file = diagn_file;
    end
elseif (nargin == 1)
    % all good - nothing to do
else
    error('Must supply exactly one filename or none.')
end

pinfo = CheckModel(input_file); % also gets default values for this model.

if (exist(truth_file,'file')==2)
   pinfo  = rmfield(pinfo,{'time_series_length','time'});
   MyInfo = CheckModelCompatibility(truth_file, input_file);
   pinfo  = CombineStructs(pinfo,MyInfo);
else
   truth_file = [];
end

pinfo.diagn_file = input_file;
pinfo.truth_file = truth_file;

clear MyInfo mynames myname ifield

%% For each model, do what needs to be done.

switch lower(pinfo.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale','lorenz_96_tracer_advection' ...
	 'forced_lorenz_96','lorenz_04','ikeda','simple_advection', 'null'}

      varid          = SetVariableID(pinfo);
      pinfo.var      = varid.var;
      pinfo.var_inds = varid.var_inds;

      fprintf('Comparing %s and \n          %s\n', pinfo.truth_file, pinfo.diagn_file)
      fprintf('Using Variable %s IDs %s\n', pinfo.var,num2str(pinfo.var_inds))
      clear varid

   case {'fms_bgrid'}

      pinfo = GetBgridInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   case {'cam'}

      pinfo.prior_file     = [];
      pinfo.posterior_file = [];
      pinfo                = GetCamInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   case {'wrf'}

      pinfo = GetWRFInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   case {'pe2lyr'}

      pinfo = GetPe2lyrInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   case {'mitgcm_ocean'}

      pinfo = GetMITgcm_oceanInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   case {'mpas_atm'}

      pinfo = GetMPAS_ATMInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   case {'sqg'}

      pinfo = GetSqgInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   case {'pop'}

      pinfo = GetPOPInfo(pinfo, input_file, 'PlotEnsTimeSeries');

   otherwise

      error('model %s not implemented yet', pinfo.model)

end

pinfo

PlotEnsTimeSeries( pinfo )


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
