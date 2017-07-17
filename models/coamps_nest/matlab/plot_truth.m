% Given a NetCDF file corresponding to the DART-generated true
% state (TrueState.nc), plot the field specified for the number of
% times specified

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear all;
close all;

% Path where the files can be found
path='/net/ds-0b/export/ensembles/tropical/3km_eakf_control';

% Determine the horizontal grid size
xdim = 100;
ydim = 100;

% Subplot size
sub_x = 3;
sub_y = 4;

% Open up the NetCDF file
filename='True_State.nc';
truth = netcdf(sprintf('%s/%s',path,filename));

times = squeeze(truth{'time'}(:));

% Only consider a single state variable - note that the first two
% items are the ensemble mean and spread
avail_times = 1;
ijarea = xdim*ydim;
total_members = 1;

% Calculate the vector position for this particular state variable
restart_info = sprintf('%s/%s',path,'restart.vars');
var_name    = 'qv2';
sigma_level = 40;
varnum = get_varnum(restart_info,var_name,sigma_level);
spots = get_state_indices(varnum,ijarea);

% Get the horizontal coordinate information
[lon,lat] = read_location(truth,spots);
lat_field = reshape(lat,[ydim xdim]);
lon_field = reshape(lon,[ydim xdim]);

% Begin figure creation
disp('Creating truth plot...')
truthfig = figure('Visible','on');


for tt=1:length(times)
  % Read in the state from the file.
  truth_state = read_state(truth,tt,1,spots);

  % Convert the vector to something we can easily plot
  % NEED TO DO THIS AS [ydim xdim] when using lat/lon!!!
  truth_state_field = reshape(truth_state,[ydim xdim]);

  % Plot the ensemble spread, mean, and truth
  subplot(sub_x,sub_y,tt)
  pcolor(lon_field*1E-3,lat_field*1E-3,truth_state_field);
  shading flat;
  axis square;
  %title(['\tau = ' sprintf('%4.2f',times(tt))]);
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
