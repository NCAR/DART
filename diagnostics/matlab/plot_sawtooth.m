function plot_sawtooth(posterior_file,prior_file)
%% DART:plot_sawtooth - time series of a state variable including updates.
%
% plot_sawtooth    interactively queries for the needed information.
%              Since different models potentially need different pieces of
%              information ... the model types are determined and additional
%              user input may be queried.
%
% Both the prior and posterior estimates are plotted as a single
% trajectory. If there is no change in the model state, this should
% appear as a series of steps. This necessitates plotting the 'posterior'
% first ... think about it ...
% If the true state ('true_state.nc') is available, it is also plotted.
%
% Example 1:  prompt for filenames 
% plot_sawtooth

% Example 2:
% posterior_file = 'analysis.nc';
% prior_file = 'preassim.nc';
% plot_sawtooth(posterior_file,prior_file)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin == 0)
    disp('Input name of (posterior) ensemble trajectory file:')
    posterior_file = input('<cr> for analysis.nc\n','s');
    if isempty(posterior_file)
        posterior_file = 'analysis.nc';
    end
    
    disp('Input name of (prior) ensemble trajectory file:')
    prior_file = input('<cr> for preassim.nc\n','s');
    if isempty(prior_file)
        prior_file = 'preassim.nc';
    end
elseif (nargin == 2)
    % proceed as normal
else
    error('Requires exactly 2 input filenames.')
end

% CheckModelCompatibility assumes first file is 'truth', so the
% components must be renamed in this context.
vars     = CheckModel(prior_file);
pinfo    = CheckModelCompatibility(prior_file, posterior_file);
pinfo.prior_time     = pinfo.truth_time;
pinfo.posterior_time = pinfo.diagn_time;
pinfo.truth_file     = 'true_state.nc';
pinfo.prior_file     = prior_file;
pinfo.posterior_file = posterior_file;
pinfo = rmfield(pinfo,{'diagn_file','truth_time','diagn_time'});
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(prior_file);
pinfo = CombineStructs(vars, pinfo);

switch lower(pinfo.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection','lorenz_96_tracer_advection','null'}

      % determine which variable ['state','X','Y' ...] and
      % determine which variable IDs (locations), as well as
      % determine which ensemble members to plot.

      pinfo                = SetVariableID(pinfo);
      pinfo.copyindices    = SetCopyID(prior_file);

   case {'fms_bgrid'}

      pinfo = GetBgridInfo(pinfo, prior_file, 'PlotSawtooth');

   case {'pe2lyr'}

      pinfo = GetPe2lyrInfo(pinfo, prior_file, 'PlotSawtooth');

   case {'wrf'}

      pinfo = GetWRFInfo(pinfo, prior_file, 'PlotSawtooth');

   case {'cam'}

      pinfo = GetCamInfo(pinfo, prior_file, 'PlotSawtooth');
      pinfo.copyindices = SetCopyID2(pinfo.prior_file);
      pinfo.copies      = length(pinfo.copyindices);

   case {'mitgcm_ocean'}

      pinfo = GetMITgcm_oceanInfo(pinfo, prior_file, 'PlotSawtooth');
      pinfo.copyindices = SetCopyID2(pinfo.prior_file);
      pinfo.copies      = length(pinfo.copyindices);

   case {'mpas_atm'}

      pinfo = GetMPAS_ATMInfo(pinfo, prior_file, 'PlotSawtooth');

   case {'sqg'}

      pinfo = GetSqgInfo(pinfo, prior_file, 'PlotSawtooth');

   case {'pop'}

      pinfo = GetPOPInfo(pinfo, prior_file, 'PlotSawtooth');

   otherwise

      error('model %s not implemented yet', pinfo.model)

end

pinfo

PlotSawtooth( pinfo )


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
