%% DART:plot_smoother_err - summary plots of global error and spread using the smoother
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
% Example 1
% diagn_file = 'postassim.nc';
% truth_file = 'perfect_output.nc';   % for smoother, is Lag_00001_Diag.nc better?
% num_lags   = 10;
% plot_total_err

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

lag_file   = 'Lag_%05d_Diag.nc'; % pattern for lag file names

if (exist('num_lags','var')   ~= 1), num_lags = 10000; end
if (exist('truth_file','var') == 1),
   def_true = truth_file;
else
   def_true = 'True_State.nc';
end

disp('Input name of True State file:')
truth_file = input(sprintf('<cr> for %s\n',def_true),'s');
if isempty(truth_file)
   truth_file = def_true;
end

if ( exist(truth_file,'file') ~= 2 ), error('%s does not exist.',truth_file); end

% Loop over all possible lags, if the corresponding etCDF file
% does not exist, we automatically terminate.

for lag=1:num_lags

  def_diag = sprintf(lag_file, lag);

  disp('Input name of smoother lag diagnostics file:')
  diagn_file = input(sprintf('<cr> for %s\n', def_diag),'s');
  if isempty(diagn_file)
     diagn_file = def_diag;
  end

  if (exist(diagn_file,'file') ~= 2)
     disp('file does not exist. Must be done.')
     return
  end

  pinfo = CheckModel(diagn_file);
  pinfo.truth_file = truth_file;
  pinfo.diagn_file = diagn_file;

  bob = CheckModelCompatibility(truth_file, diagn_file);
  pinfo.truth_time = bob.truth_time;
  pinfo.diagn_time = bob.diagn_time;

  clear bob

  fprintf('Comparing %s and \n          %s\n', ...
                pinfo.truth_file, pinfo.diagn_file)

  PlotTotalErr( pinfo );

  disp(' ')

end

clear pinfo


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
