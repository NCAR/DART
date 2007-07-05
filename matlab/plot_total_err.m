% DART : summary plots of global error and spread
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Posterior_Diag.nc';
% plot_total_err

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
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

pinfo = struct('truth_file', truth_file, ...
               'diagn_file', diagn_file, ...
               'truth_time', [-1,-1], ...
               'diagn_time', [-1,-1] );

disp(sprintf('Comparing %s and \n          %s', truth_file, diagn_file))

PlotTotalErr( pinfo )
clear pinfo
