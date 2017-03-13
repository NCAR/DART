function plot_ens_err_spread(truth_file, diagn_file)
%% DART: plot_ens_err_spread - summary plots of the ensemble error and ensemble spread.
%                        Interactively queries for the needed information.
%                        Since different models potentially need different
%                        pieces of information ... the model types are
%                        determined and additional user input may be queried.
%
% Ultimately, plot_ens_err_spread will be replaced by a GUI.
% All the heavy lifting is done by PlotEnsErrSpread.
%
% A reminder of the sequence:
% truth  run (from    pmo):
%           perfect_input  --->  perfect_output.nc
% filter run (from filter):
%           filter_input.nc  --->  [prior inflation]  --->
%                 preassim.nc   --->  [assimilation]  --->
%                       postassim.nc  ---> [posterior inflation]  --->
%                             filter_output.nc
%
% Example 1  (Prompt for filenames. Defaults are 'perfect_output.nc' and 'preassim.nc')
% plot_ens_err_spread
%
% Example 2
% truth_file = 'perfect_output.nc';
% diagn_file = 'preassim.nc';
% plot_ens_err_spread(truth_file, diagn_file)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin == 0)
    disp('Input name of true model trajectory file:')
    truth_file = input('<cr> for perfect_output.nc\n','s');
    if isempty(truth_file)
        truth_file = 'perfect_output.nc';
    end
    disp('Input name of ensemble trajectory file:')
    diagn_file = input('<cr> for preassim.nc\n','s');
    if isempty(diagn_file)
        diagn_file = 'preassim.nc';
    end
elseif (nargin == 2)
    % all good - nothing to do
else
    error('Must supply either two filenames or none.')
end

if ( exist(truth_file,'file') ~= 2 ), error('%s does not exist.',truth_file); end
if ( exist(diagn_file,'file') ~= 2 ), error('%s does not exist.',diagn_file); end

vars  = CheckModel(diagn_file);   % also gets default values for this model.
vars  = rmfield(vars,{'fname','time','time_series_length'});
pinfo = CheckModelCompatibility(truth_file, diagn_file);
pinfo = CombineStructs(pinfo,vars);
clear vars

switch lower(pinfo.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection'}

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
