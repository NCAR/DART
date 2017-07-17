% Plot 2-dimensional field maps of the ensemble spread for either
% the prior or posterior diagnostic NetCDF file.  

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear all;
close all;

% Path where the files can be found
path='/net/ds-0b/export/ensembles/tropical/30km_eakf_control';

% Determine the horizontal grid size
xdim = 100;
ydim = 100;

% Open up the NetCDF file
filename='Prior_Diag.nc';
diag_nc = netcdf(sprintf('%s/%s',path,filename));

% Grab the entire set of times
times = squeeze(diag_nc{'time'}(:));

% Calculate the position in the big state vector for this
% particular variable 
restart_info = sprintf('%s/%s',path,'restart.vars');
var_name    = 'v2';
sigma_level = 40;
varnum = get_varnum(restart_info,var_name,sigma_level);
ijarea = xdim*ydim;
spots = get_state_indices(varnum,ijarea);

% Get the horizontal coordinate information
[lon,lat] = read_location(diag_nc,spots);
lat_field = reshape(lat,[xdim ydim]).';
lon_field = reshape(lon,[xdim ydim]).';

% Generate figures offline
disp('Creating spread plot...')
diag_ncfig = figure('Visible','off');

% Generate one figure for each time
for tt=1:length(times)
  % Read in the spread - this is element 2 
  spread = read_state(diag_nc,tt,2,spots);

  % Convert the vector to something we can easily plot
  spread_field = reshape(spread,[xdim ydim]).';

  % Generate the color plot
  pcolor(lon_field,lat_field,spread_field);
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

saveas(diag_ncfig,'spread','png');

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
