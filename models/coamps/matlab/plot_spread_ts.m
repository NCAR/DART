%% plot_spread_ts   plot the spread field specified for the number of times specified.
%
% Uses a DART-generated NetCDF file: true_state.nc, preassim.nc, analysis.nc

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear all;
close all;

% Path where the files can be found
path='/net/ds-0b/export/ensembles/tropical/3km_control';

% Determine the horizontal grid size
xdim = 100;
ydim = 100;

% Open up the NetCDF file
filename='analysis.nc';
diag_nc = netcdf(sprintf('%s/%s',path,filename));

% Grab the entire set of times
times = squeeze(diag_nc{'time'}(:));

% Calculate the position in the big state vector for this
% particular variable 
restart_info = sprintf('%s/%s',path,'restart.vars');
var_name    = 'qv2';
sigma_level = 45;
varnum = get_varnum(restart_info,var_name,sigma_level);
ijarea = xdim*ydim;
spots = get_state_indices(varnum,ijarea);

% Grab the spread at each time
for tt=1:length(times)
  spread = read_state(diag_nc,tt,2,spots);
  total_spread(tt) = sum(spread);
end

% Begin figure creation - note that time is in fractional days
disp('Creating spread plot...')
diag_ncfig = figure('Visible','off');
plot(times*24,total_spread);
xlabel('Time (h)');
ylabel('Ensemble Spread');
title([var_name ',\sigma=' sprintf('%d',sigma_level)]);
saveas(diag_ncfig,'spread','png');

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
