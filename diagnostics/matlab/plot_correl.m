function plot_correl(diagn_file)
%% DART:plot_correl Plots space-time series of correlation between a given variable
%               at a given time and other variables at all times in an
%               ensemble time sequence.
%
% plot_correl  interactively queries for the information needed to create
%              the desired correlations.
%              Since different models potentially need different pieces
%              of information ... the model types are determined and
%              additional user input may be queried.
%
% A reminder of the sequence:
% truth  run (from    pmo):
%           perfect_input  --->  true_state.nc
% filter run (from filter):
%           filter_input.nc  --->  [prior inflation]  --->
%                 preassim.nc   --->  [assimilation]  --->
%                       postassim.nc  ---> [posterior inflation]  --->
%                             filter_output.nc
%
% Example 1 (prompts you for a filename - default is 'preassim.nc')
% plot_correl
%
% Example 2
% fname = 'analysis.nc';
% plot_correl(fname)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin == 0)
    disp('Input name of ensemble trajectory file:')
    diagn_file = input('<cr> for preassim.nc\n','s');
    if isempty(diagn_file)
        diagn_file = 'preassim.nc';
    end
elseif (nargin ==1 )
    % nothing to do
else
    error('Must supply 1 filename or none.')
end

if ( exist(diagn_file,'file') ~= 2 ), error('%s does not exist.',diagn_file); end

pinfo = CheckModel(diagn_file);

% check to make sure they are using a file with some ensemble members.
if (pinfo.num_ens_members <= 3)
   error('Sorry -- %s does not have enough ensemble members to correlate.',diagn_file)
end

switch lower(pinfo.model)
   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04', ...
	 'forced_lorenz_96','ikeda', 'null'}

      pinfo.base_var = pinfo.def_var;

      inputstring = input( ...
           sprintf('Input index for base variable (between %d and %d, inclusive)  ', ...
           pinfo.min_state_var,pinfo.max_state_var),'s');
      pinfo.base_var_index = str2double(deblank(inputstring));

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d, inclusive)  ', ...
           pinfo.time_series_length),'s');
      pinfo.base_time = str2double(deblank(inputstring));

      fprintf('Using diagnostic file %s\n',diagn_file)
      fprintf('Correlating variable %s index %d at time %d.\n', ...
           pinfo.base_var, pinfo.base_var_index, pinfo.base_time)

   case {'lorenz_96_2scale'}

      fprintf('Your choice of variables is ''X'' or ''Y''\n')
      fprintf('''X'' can range from %d to %d\n', pinfo.min_X_var, pinfo.max_X_var)
      fprintf('''Y'' can range from %d to %d\n', pinfo.min_Y_var, pinfo.max_Y_var)

      % parsing the result of this one is a bit tricky.
      inputstring = input('Input base variable and index i.e.  X 5\n','s');
      [pinfo.base_var, pinfo.base_var_index] = ParseAlphaNumerics(inputstring);

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d, inclusive)  ', ...
           pinfo.time_series_length),'s');
      pinfo.base_time = str2double(deblank(inputstring));

      fprintf('Using diagnostic file %s\n',diagn_file)
      fprintf('Correlating variable %s index %d at time %d.\n', ...
           pinfo.base_var,pinfo.base_var_index, pinfo.base_time)

   case {'simple_advection', 'lorenz_96_tracer_advection'}

      disp('Your choice of variables are:')
      disp(pinfo.vars)
      fprintf('the indices (locations) can range from %d to %d, inclusive\n', ...
           pinfo.min_state_var, pinfo.max_state_var)

      str1 = sprintf('Input base variable and index i.e. %s %d\n', ...
                      pinfo.def_var,pinfo.def_state_vars(1));
      inputstring = input(str1,'s');
      [pinfo.base_var, pinfo.base_var_index] = ParseAlphaNumerics(inputstring);

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d, inclusive)  ', ...
           pinfo.time_series_length),'s');
      pinfo.base_time = str2double(deblank(inputstring));

      fprintf('Using diagnostic file %s\n',diagn_file)
      fprintf('Correlating variable %s index %d at time %d.\n', ...
           pinfo.base_var,pinfo.base_var_index, pinfo.base_time)

   case {'fms_bgrid'}

      pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'cam'}

      pinfo = GetCamInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'clm'}

      pinfo = GetClmInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'wrf'}

      pinfo = GetWRFInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'pe2lyr'}

      pinfo = GetPe2lyrInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'mitgcm_ocean'}

      pinfo = GetMITgcm_oceanInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'tiegcm'}

      pinfo = GetTIEGCMInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'mpas_atm'}

      pinfo = GetMPAS_ATMInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'sqg'}

      pinfo = GetSqgInfo(pinfo, diagn_file, 'PlotCorrel');

   case {'pop'}

      pinfo = GetPOPInfo(pinfo, diagn_file, 'PlotCorrel');

   otherwise

      error('model %s not implemented yet', pinfo.model)

end

pinfo

PlotCorrel( pinfo );
clear inputstring inds str1 vrbl vrbl_inds diminfo


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
