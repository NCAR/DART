%% plot_coamps_nogfx -  Makes some diagnostic plots for DART conversion
% ----------
% plot_l05.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear all;
close all;

% Open up the various NetCDF files
truth = netcdf('true_state.nc');
prior = netcdf('preassim.nc');
postr = netcdf('analysis.nc');

% Pull out the time variable
state_time = prior{'time'}(:);

% Pull out the state vars in each one
truth_fullstate = truth{'state'}(:);
prior_fullstate = prior{'state'}(:);
postr_fullstate = postr{'state'}(:);

total_members = 10;

% Only consider a single state variable - note that the first two
% items are the ensemble mean and spread
avail_times = 2;
member = 6;
ijarea = 72*45;
varnum = 2;
spots = ((varnum-1)*ijarea+1):(varnum*ijarea);
truth_state  = squeeze(truth_fullstate(avail_times,:,spots));
prior_spread = squeeze(prior_fullstate(avail_times,2,spots));
prior_mean   = squeeze(prior_fullstate(avail_times,1,spots));
postr_spread = squeeze(postr_fullstate(avail_times,2,spots));
postr_mean   = squeeze(prior_fullstate(avail_times,1,spots));

for ii=[1:total_members]
  prior_state(:,:,ii) = squeeze(prior_fullstate(avail_times,ii,spots));
  postr_state(:,:,ii) = squeeze(postr_fullstate(avail_times,ii,spots));
end

% Convert the vector to something we can easily plot
truth_state_field = reshape(truth_state,[72 45]).';

for ii=1:total_members
  prior_state_field(:,:,ii) = reshape(prior_state(:,:,ii),[72 45]).';
  postr_state_field(:,:,ii) = reshape(postr_state(:,:,ii),[72 45]).';
end

prior_spread_field = reshape(prior_spread,[72 45]).';
postr_spread_field = reshape(postr_spread,[72 45]).';
prior_mean_field   = reshape(  prior_mean,[72 45]).';
postr_mean_field   = reshape(  postr_mean,[72 45]).';

% standard contour levels
stdlevels = linspace(min(truth_state),max(truth_state),10);

% Plot the ensemble spread, mean, and truth
priorfig = figure('Visible','off');
contourf(prior_spread_field);
cur_ax = gca;
caxis(cur_ax,[min(min(prior_spread_field)) max(max(prior_spread_field))]);
shading flat;
colorbar('h');
hold on;
for ii=1:total_members
  [c_mem,h_mem] = contour(prior_state_field(:,:,ii),stdlevels);
  set(h_mem, 'LineWidth',1.0,'LineStyle','-','LineColor',[1 0 0])
end
[c_truth,h_truth] = contour(truth_state_field,stdlevels, 'LineColor',...
                            [0 0 0]);
axis image;
set(h_truth,'LineWidth',1.0,'LineStyle','-','LineColor',[0 1 0]);
xlabel('X Gridpoint');
ylabel('Y Gridpoint');
title (['Prior Estimate of Ensemble MSLP Mean and Spread: 31 August 1999 : \tau = ' sprintf('%dhrs',avail_times ...
                                           + 5)]);
leg = legend([h_mem(1) h_truth(1)],'Ensemble Members','Truth', ...
             'Location','SouthOutside');
saveas(priorfig,'prior','png');

% Plot the ensemble spread, mean, and truth
postrfig = figure('Visible','off');
contourf(postr_spread_field);
cur_ax = gca;
caxis(cur_ax,[min(min(postr_spread_field)) max(max(postr_spread_field))]);
shading flat;
colorbar('h');
hold on;
for ii=1:total_members
  [c_mem,h_mem] = contour(postr_state_field(:,:,ii),stdlevels);
  set(h_mem, 'LineWidth',1.0,'LineStyle','-','LineColor',[1 0 0])
end
[c_truth,h_truth] = contour(truth_state_field,stdlevels, 'LineColor',...
                            [0 0 0]);
axis image;
set(h_truth,'LineWidth',1.0,'LineStyle','-','LineColor',[0 1 0]);
xlabel('X Gridpoint');
ylabel('Y Gridpoint');
title (['Posterior Estimate of Ensemble MSLP Mean and Spread: 31 August 1999 : \tau = ' sprintf('%dhrs',avail_times ...
                                           + 5)]);
leg = legend([h_mem(1) h_truth(1)],'Ensemble Members','Truth', ...
             'Location','SouthOutside');
saveas(postrfig,'posterior','png');

% Plot the ensemble spread, mean, and truth
difffig = figure('visible','off');
diff_field = (postr_state_field(:,:,member)-prior_state_field(:,:,member));
contourf(diff_field);
shading flat;
colorbar('h')
cur_ax = gca;
hold on;
axis image;
xlabel('X Gridpoint');
ylabel('Y Gridpoint');
title (['Difference between Prior and Posterior Estimates: ' ...
        '31 August 1999 : \tau = ' sprintf('%dhrs',avail_times ...
                                           + 5)]);
saveas(difffig,'difference','png');

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
