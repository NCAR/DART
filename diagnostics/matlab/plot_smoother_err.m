function plot_smoother_err
%% DART:plot_smoother_err - summary plots of global error and spread using the smoother
%
% Example 1
% diagn_file = 'postassim.nc';
% num_lags   = 10;
% plot_total_err

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

lag_file   = 'Lag_%05d_Diag.nc'; % pattern for lag file names

if (exist('num_lags','var')   ~= 1), num_lags = 10000; end
truth_file = 'true_state.nc';

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
