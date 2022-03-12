function plot_ens_err_spread(input_file)
%% DART: plot_ens_err_spread - summary plots of the ensemble error and ensemble spread.
%                        Interactively queries for the needed information.
%                        Since different models potentially need different
%                        pieces of information ... the model types are
%                        determined and additional user input may be queried.
%
% The true state file is REQUIRED and is specified by the variable 'truth_file'.
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
% plot_ens_err_spread
%
% Example 2:  Prompt, but use name specified in 'diagn_file'.
% diagn_file = 'analysis.nc';
% plot_ens_err_spread
%
% Example 3:  Prompt with 'diagn_file', different true state.
% diagn_file = 'analysis.nc';
% truth_file = 'More_True_State.nc';
% plot_ens_err_spread
%
% Example 4: No prompting, true state filename comes from 'truth_file' variable.
% plot_ens_err_spread('output.nc')

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

if ( exist(truth_file,'file') ~= 2 ), error('%s does not exist.',truth_file); end
if ( exist(input_file,'file') ~= 2 ), error('%s does not exist.',input_file); end

vars  = CheckModel(input_file);   % also gets default values for this model.
vars  = rmfield(vars,{'fname','time','time_series_length'});
pinfo = CheckModelCompatibility(truth_file, input_file);
pinfo = CombineStructs(pinfo,vars);
clear vars

%% For each model, do what needs to be done.

switch lower(pinfo.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'forced_lorenz_96','lorenz_04','ikeda','simple_advection','lorenz_96_tracer_advection', 'null'}

      varid          = SetVariableID(pinfo);
      pinfo.var      = varid.var;
      pinfo.var_inds = varid.var_inds;

      fprintf('Comparing %s and \n          %s\n', pinfo.truth_file, pinfo.diagn_file)
      fprintf('Using Variable %s IDs %s\n', pinfo.var,num2str(pinfo.var_inds))
      clear varid

   case {'fms_bgrid'}

      pinfo = GetBgridInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case {'cam'}

      pinfo = GetCamInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case {'wrf'}

      pinfo = GetWRFInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case {'pe2lyr'}

      pinfo = GetPe2lyrInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case {'mitgcm_ocean'}

      pinfo = GetMITgcm_oceanInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case {'mpas_atm'}

      pinfo = GetMPAS_ATMInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case {'sqg'}

      pinfo = GetSqgInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   case {'pop'}

      pinfo = GetPOPInfo(pinfo, truth_file, 'PlotEnsErrSpread');

   otherwise

      error('model %s not implemented yet', pinfo.model)

end

pinfo

PlotEnsErrSpread( pinfo )


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
