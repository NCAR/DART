% DART : summary plots of global error and spread using the smoother
% Example 1
% diagn_file = 'Posterior_Diag.nc';
% truth_file = 'True_State.nc';   % for smoother, is Lag_00001_Diag.nc better?
% num_lags   = 10;
% plot_total_err

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL: http://subversion.ucar.edu/DAReS/DART/trunk/matlab/plot_total_err.m $
% $Id: plot_total_err.m 2691 2007-03-11 18:18:09Z thoar $
% $Revision: 2691 $
% $Date: 2007-03-11 12:18:09 -0600 (Sun, 11 Mar 2007) $

% pattern for how lag files are written out
lag_file   = 'Lag_%05d_Diag.nc';

num_lags = 10000

for lag=1:num_lags
  
  if (exist(truth_file) ~= 2)
     disp('Input name of Truth State file;')
     truth_file = input('<cr> for True_State.nc\n','s');
     if isempty(truth_file)
        truth_file = 'True_State.nc';
     end
  end

  diagn_file = sprintf(lag_file, lag)
  prompt = sprintf('<cr> for %s\n', diagn_file);
  
  if (exist(diagn_file) ~=2)
     disp('Input name of smoother lag diagnostics file;')
     diagn_file = input(prompt,'s');
     if isempty(diagn_file)
        return    % nothing more to do
     end
  end
  
  
  pinfo = CheckModelCompatibility(truth_file, diagn_file)
  
  disp(sprintf('Comparing %s and \n          %s', ...
                pinfo.truth_file, pinfo.diagn_file))
  
  %disp('in main after return')
  %pinfo
  
  PlotTotalErr( pinfo );

  input('Hit return to continue');
end

clear pinfo
