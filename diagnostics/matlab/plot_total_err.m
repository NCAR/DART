function plot_total_err(input_file)
%% DART:plot_total_err - summary plots of global error and spread
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
% plot_total_err
%
% Example 2:  Prompt, but use name specified in 'diagn_file'.
% diagn_file = 'analysis.nc';
% plot_total_err
%
% Example 3:  Prompt with 'diagn_file', different true state.
% diagn_file = 'analysis.nc';
% truth_file = 'More_True_State.nc';
% plot_total_err
%
% Example 4: No prompting, true state filename comes from 'truth_file' variable.
% plot_total_err('output.nc')

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

fprintf('Comparing %s and \n          %s\n', truth_file, input_file)

CheckModel(input_file); % make sure model is supported - no need for anything else.
pinfo   = CheckModelCompatibility(truth_file,input_file);

switch lower(pinfo.model)
    case{'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale','lorenz_96_tracer_advection', ...
         'lorenz_04', 'forced_lorenz_96','ikeda','simple_advection', 'null'}

    case{'fms_bgrid'}
        pinfo = GetBgridInfo(pinfo, input_file, 'PlotTotalErr');

    case {'wrf'}
      pinfo = GetWRFInfo(pinfo, input_file, 'PlotTotalErr');

    case {'cam'}
      pinfo = GetCamInfo(pinfo, input_file, 'PlotTotalErr');

   case {'pe2lyr'}
      pinfo = GetPe2lyrInfo(pinfo, input_file, 'PlotTotalErr');

   case {'mitgcm_ocean'}
      pinfo = GetMITgcm_oceanInfo(pinfo, input_file, 'PlotTotalErr');

   case {'mpas_atm'}
      pinfo = GetMPAS_ATMInfo(pinfo, input_file, 'PlotTotalErr');

   case {'sqg'}
      pinfo = GetSqgInfo(pinfo, input_file, 'PlotTotalErr');

   case {'pop'}
      pinfo = GetPOPInfo(pinfo, input_file, 'PlotTotalErr');

   otherwise

      error('%s not implemented yet', pinfo.model)
end

PlotTotalErr( pinfo );


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
