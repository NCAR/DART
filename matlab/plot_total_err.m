%% DART:plot_total_err - summary plots of global error and spread
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Posterior_Diag.nc';
% plot_total_err

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist('truth_file','var') ~= 1)
   disp('Input name of True State file:')
   truth_file = input('<cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('diagn_file','var') ~=1)
   disp('Input name of prior or posterior diagnostics file:')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

if ( exist(truth_file,'file') ~= 2 ), error('%s does not exist.',truth_file); end
if ( exist(diagn_file,'file') ~= 2 ), error('%s does not exist.',diagn_file); end

fprintf('Comparing %s and \n          %s\n', truth_file, diagn_file)

CheckModel(diagn_file); % make sure model is supported - no need for anything else.
pinfo   = CheckModelCompatibility(truth_file,diagn_file);

switch lower(pinfo.model)
    case{'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
         'lorenz_04', 'forced_lorenz_96','ikeda','simple_advection'}

    case{'fms_bgrid'}
        pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotTotalErr');

    case {'wrf'}
      pinfo = GetWRFInfo(pinfo, diagn_file, 'PlotTotalErr');

    case {'cam'}
      pinfo = GetCamInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'pe2lyr'}
      pinfo = GetPe2lyrInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'mitgcm_ocean'}
      pinfo = GetMITgcm_oceanInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'mpas_atm'}
      pinfo = GetMPAS_ATMInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'sqg'}
      pinfo = GetSqgInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'pop'}
      pinfo = GetPOPInfo(pinfo, diagn_file, 'PlotTotalErr');

   otherwise

      error('%s not implemented yet', pinfo.model)
end

PlotTotalErr( pinfo );
clear pinfo


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

