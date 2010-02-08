%% DART:plot_total_err - summary plots of global error and spread
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Posterior_Diag.nc';
% plot_total_err

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist('truth_file','var') ~= 1)
   truth_file = input('Input name of True State file; <cr> for True_State.nc\n','s');
   if isempty(truth_file)
      truth_file = 'True_State.nc';
   end
end

if (exist('diagn_file','var') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end

if ( exist(truth_file,'file') ~= 2 ), error('%s does not exist.',truth_file); end
if ( exist(diagn_file,'file') ~= 2 ), error('%s does not exist.',diagn_file); end

disp(sprintf('Comparing %s and \n          %s', truth_file, diagn_file))

pinfo = CheckModel(diagn_file);
pinfo.truth_file = truth_file;
pinfo.diagn_file = diagn_file;
pinfo.truth_time = [1, -1];    % The temporal intersection happens in PlotTotalErr.m
pinfo.diagn_time = [1, -1];

PlotTotalErr( pinfo );
clear pinfo
