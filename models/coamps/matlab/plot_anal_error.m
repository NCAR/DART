%% Plot 2-dimensional field maps of the analysis error for either
% the prior or posterior diagnostic NetCDF file.  

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear all;
close all;

% Path where the files can be found
path='/net/ds-01/scratch-5/timw/nrl/unconstrained';

% Determine the horizontal grid size
xdim = 92;
ydim = 55;

% Open up the NetCDF files - need the ensemble mean and truth
diag_filename='analysis.nc';
diag_nc = netcdf(sprintf('%s/%s',path,diag_filename));
true_filename='true_state.nc';
true_nc = netcdf(sprintf('%s/%s',path,true_filename));

% Grab the entire set of times
times = squeeze(true_nc{'time'}(:));

% Calculate the position in the big state vector for this
% particular variable 
restart_info = sprintf('%s/%s',path,'restart.vars');
var_name    = 'u2';
sigma_level = 25;
varnum = get_varnum(restart_info,var_name,sigma_level);
ijarea = xdim*ydim;
spots = get_state_indices(varnum,ijarea);

% Get the horizontal coordinate information
[lon,lat] = read_location(diag_nc,spots);
lat_field = reshape(lat,[xdim ydim]).';
lon_field = reshape(lon,[xdim ydim]).';

% Generate figures offline
disp('Creating analysis error plot...')
analerr_fig = figure('Visible','off');

% Generate one figure for each time
for tt=1:length(times)
  % Read in the true state
  truth = read_state(true_nc,tt,1,spots);
  
  % Read in the ensemble mean as the best estimate
  estimate = read_state(diag_nc,tt,1,spots);

  % We want the magnitude of the error: Since we're plotting the
  % timeseries as RMSE, do the same style of computation here.
  error_mag = sqrt((truth - estimate).^2);
  
  % Convert the vector to something we can easily plot
  error_field = reshape(error_mag,[xdim ydim]).';

  % Generate the color plot
  pcolor(lon_field,lat_field,error_field);
  shading flat;
  axis image;
  
  % These labels are straightforward
  xlabel('Longitude ( ^o E)');
  ylabel('Latitude ( ^o N)');
  
  % Label with the variable name, the sigma level, and the forecast
  % time.  In DART, forecast time is in days, so convert to hours 
  title([var_name ',\sigma=' sprintf('%d',sigma_level) ...
         '\tau = ' sprintf('%4.2f',times(tt)*24)]);
end

saveas(analerr_fig,sprintf('anal_error_%s_%d',var_name,sigma_level),'png');

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
