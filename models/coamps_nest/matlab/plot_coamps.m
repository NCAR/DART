%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear all;
close all;


xdim = 28;
ydim = 28;

% Open up the various NetCDF files
%truth = netcdf('true_state.nc');
prior = netcdf('preassim.nc');
postr = netcdf('analysis.nc');


% Only consider a single state variable - note that the first two
% items are the ensemble mean and spread
avail_times = 3;
ijarea = xdim*ydim;
varnum = 372;
total_members = 100;

% Calculate the vector position for this particular state variable field
spots = ((varnum-1)*ijarea+1):(varnum*ijarea);

% Read in the true value, the prior and posterior mean/spread
%truth_state = squeeze(truth{'state'}(avail_times,:,spots));
prior_spread = squeeze(prior{'state'}(avail_times,2,spots));
prior_mean = squeeze(prior{'state'}(avail_times,1,spots));
postr_spread = squeeze(postr{'state'}(avail_times,2,spots));
postr_mean = squeeze(prior{'state'}(avail_times,1,spots));

% Have to iterate to get all the ensemble members - offset by 2
% since the mean and spread take up the first slpts
for ii=[1:total_members]
  prior_state(:,:,ii) = squeeze(prior{'state'}(avail_times,ii+2,spots));
  postr_state(:,:,ii) = squeeze(postr{'state'}(avail_times,ii+2,spots));
end

% Lat/Lon information comes from the truth run
lat = squeeze(prior{'lat'}(spots));
lon = squeeze(prior{'lon'}(spots));

% Convert the vector to something we can easily plot
% NEED TO DO THIS AS [ydim xdim] when using lat/lon!!!
%truth_state_field = reshape(truth_state,[ydim xdim]);
for ii=1:total_members
  prior_state_field(:,:,ii) = reshape(prior_state(:,:,ii),[ydim xdim]);
  postr_state_field(:,:,ii) = reshape(postr_state(:,:,ii),[ydim xdim]);
end
prior_spread_field = reshape(prior_spread,[ydim xdim]);
postr_spread_field = reshape(postr_spread,[ydim xdim]);
prior_mean_field = reshape(prior_mean,[ydim xdim]);
postr_mean_field = reshape(postr_mean,[ydim xdim]);

lat_field = reshape(lat,[ydim xdim]);
lon_field = reshape(lon,[ydim xdim]);

% Observation information
obs_spacing = 4;

lat_obs = lat_field(1:obs_spacing:end,1:obs_spacing:end);
lon_obs = lon_field(1:obs_spacing:end,1:obs_spacing:end);

% Standard contour levels
n_contours = 5;
%stdlevels =
%linspace(min(truth_state),max(truth_state),n_contours);
stdlevels = linspace(min(prior_mean),max(prior_mean),n_contours);


% Plot the ensemble spread, mean, and truth
disp('Creating prior plot...')
priorfig = figure('Visible','on');
contourf(lon_field,lat_field,prior_spread_field);
cur_ax = gca;
%caxis(cur_ax,[min(min(prior_spread_field)) max(max(prior_spread_field))]);
shading flat;
alpha(0.3);
colorbar('h');
hold on;
% for ii=1:10:total_members
%   [c_mem,h_mem] = contour(lon_field, lat_field, prior_state_field(: ...
% 						  ,:,ii),stdlevels); 
%   set(h_mem, 'LineWidth',1.0,'LineStyle','-','LineColor',[1 0 0])
% end
[c_truth,h_truth] = contour(lon_field, lat_field, ...
			    prior_mean_field,stdlevels, 'LineColor', ...
			    [0 0 0]); 
clabel(c_truth,h_truth);
axis image;
set(h_truth,'LineWidth',2.0,'LineStyle','-','LineColor',[0 0 0]);
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees N)');
%set(cur_ax,'XAxisLocation','top')
%saveas(priorfig,'prior','png');


% Plot the ensemble spread, mean, and truth
disp('Creating posterior plot...')
postrfig = figure('Visible','on');
contourf(lon_field,lat_field,postr_spread_field);
cur_ax = gca;
caxis(cur_ax,[min(min(postr_spread_field)) max(max(postr_spread_field))]);
shading flat;
%alpha(0.3);
colorbar('h');
hold on;
[c_truth,h_truth] = contour(lon_field, lat_field, ...
			    postr_mean_field,stdlevels, 'LineColor', ...
			    [0 0 0]);
clabel(c_truth,h_truth);
axis image;
set(h_truth,'LineWidth',2.0,'LineStyle','-','LineColor',[0 0 0]);
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees N)');


% Plot the difference between prior and posterior
disp('Creating difference plot...')
difffig = figure('visible','on');

disp('Calculating difference field')
diff_field = sqrt(sum((postr_state_field-prior_state_field).^2,3))/100;
disp('Plotting difference field')
pcolor(lon_field,lat_field,diff_field);
cur_ax = gca;
disp('Trying to shade...')
shading flat;
colorbar('h')
hold on;
axis image;
plot(lon_obs,lat_obs,'xk');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees N)');
%set(cur_ax,'XAxisLocation','top')
%saveas(difffig,'difference','png');


print -f2 -dtiff postr20
print -f3 -dtiff diffr20

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
