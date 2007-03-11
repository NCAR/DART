% DART : Plots time series of a state variable including updates.
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

% Ultimately, plot_sawtooth will be replaced by a GUI.
%
%
% All the heavy lifting is done by PlotSawtooth.

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

if (exist('truth_file') ~= 1)
   disp('If the True_State.nc exists, it will be plotted. If not, don''t worry.')
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('posterior_file') ~=1)
   disp('Input name of posterior diagnostics file;')
   posterior_file = input('<cr> for Posterior_Diag.nc\n','s');
   if isempty(posterior_file)
      posterior_file = 'Posterior_Diag.nc';
   end
end

if (exist('prior_file') ~=1)
   disp('Input name of prior diagnostics file;')
   prior_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(prior_file)
      prior_file = 'Prior_Diag.nc';
   end
end

CheckModelCompatibility(prior_file, posterior_file)
pstruct                = CheckModel(posterior_file);   % also gets default values
pstruct.prior_file     = prior_file;
pstruct.posterior_file = posterior_file;

if ( exist(truth_file) == 2 )
   CheckModelCompatibility(truth_file, prior_file)
   pstruct.truth_file  = truth_file;
else
   pstruct.truth_file  = [];
end


switch lower(pstruct.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96'} 

      % determine which variable ['state','X','Y' ...] and
      % determine which variable IDs (locations), as well as
      % determine which ensemble members to plot.

      pinfo               = SetVariableID(pstruct);
      pstruct.var         = pinfo.var;
      pstruct.var_inds    = pinfo.var_inds;
      pstruct.copyindices = SetCopyID(prior_file);

   case 'fms_bgrid'

      pstruct = GetBgridInfo(prior_file, 'PlotSawtooth');
      pstruct.prior_file     = prior_file;
      pstruct.posterior_file = posterior_file;

   case 'cam'

      pstruct = GetCamInfo(pstruct,'PlotSawtooth');
      pstruct.copyindices = SetCopyID2(pstruct.prior_file);
      pstruct.copies      = length(pstruct.copyindices);

   otherwise

      error(sprintf('model %s not implemented yet', pstruct.model))

end

pstruct

PlotSawtooth( pstruct )
clear pstruct pinfo
