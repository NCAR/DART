function [prior, poste] = AMS_AGU_poster
%
%
%

prior = [];
poste = [];

% H0 files have FSNO(time,lat,lon) and H2OSNO(time,lat,lon)

% /glade/scratch/thoar/archive/clm_freerun/lnd/hist 
%    clm_freerun.clm2_0031.h0.2001-01-01-00000.nc   has 31 timesteps in it ...
% mcdate ... has the (integer) YYYMMDD timestamps

[prior, poste] = Create_figure_1();
% Create_figure_3();



function [prior, poste] = Create_figure_1

fdir = '/glade/scratch/thoar/archive/clm_Tb_osse/dart/hist';

%% Figure 1 
%  plot of observation locations ... superimposed on plot of innovations.
%  difference of free run and assimilated run at 01-31-00000
%  only impacting CLM 'frac_sno' variable 'T_SOISNO'
% 'frac_sno',    'KIND_SNOWCOVER_FRAC',
% 'H20SNO',      'KIND_SNOW_WATER',
% 'H2OSOI_LIQ',  'KIND_SOIL_MOISTURE',
% 'H2OSOI_ICE',  'KIND_ICE',
% 'T_SOISNO',    'KIND_SOIL_TEMPERATURE',

figure(1); clf;
h = subplot(1,1,1);
set(h,'FontName', 'Arial', 'Fontsize', 20, 'FontWeight', 'bold')

varname    = 'frac_sno';
levelindex = 1;
timeindex  = 1;
copystring = 'ensemble mean';
clmfname   = '/glade/scratch/thoar/enstest_0907/enstest_0907.clm2_0003.r.2000-01-05-00000.nc';

priorfname = sprintf('%s/clm_Prior_Diag.2001-01-31-00000.nc',fdir);
postefname = sprintf('%s/clm_Posterior_Diag.2001-01-31-00000.nc',fdir);

prior = clm_get_var(priorfname,varname,copystring,levelindex,timeindex);
poste = clm_get_var(postefname,varname,copystring,levelindex,timeindex);

loninds = find( (prior.lon >= 210) & (prior.lon <= 310));
latinds = find( (prior.lat >=  25) & (prior.lat <=  60));
lons    = prior.lon(loninds);
lats    = prior.lat(latinds);

datmat  = poste.datmat - prior.datmat;
inds = find(datmat < -0.049);
datmat(inds) = -0.040;

subset  = datmat(min(latinds):max(latinds),min(loninds):max(loninds));
A2      = imresize(subset, 100, 'nearest');
h1      = imagesc(lons,lats,A2);
%set(h1,'AlphaData',~isnan(subset))

axis image;

set(gca,'YDir','normal')
h2 = title('Innovations in Snow Cover Fraction', ...
           'FontName','Arial','FontSize',20,'FontWeight','bold')
h3 = worldmap;
h4 = colorbar;
set(h4,'YTickLabel',['-0.05';' 0.00';' 0.05'])
set(h4,'FontName','Arial','FontSize',20,'FontWeight','bold')
ylabel('Latitude');
xlabel('Longitude');

set(gca,'Clim',[-0.05, 0.05]);
bob = jet(65);
bob(33,:) = 0.8;
bob(1,:) = 1.0;
colormap(bob)

% superimpose some observation locations ...

fname         = 'obs_epoch_001.nc';
ObsTypeString = 'AMSRE_BRIGHTNESS_T';
region        = [210 310 25 60 -Inf Inf];
CopyString    = 'observations';
QCString      = 'DART quality control';
maxgoodQC     = 2;
verbose       = 1;   % anything > 0 == 'true'
 
obs = read_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, verbose);

hold on;
h = plot(obs.lons,obs.lats,'o');
set(h,'MarkerSize',10,'LineWidth',1.5)
hold off;

set(gca,'FontName','Arial','FontSize',20,'FontWeight','bold') 
title('Innovations in Snow Cover Fraction','FontName','Arial','FontSize',20,'FontWeight','bold')


function Create_figure_2

%% Figure 2 
%  plot of trajectories 
%  difference of free run and assimilated run as a function of time.
%  what variables


function Create_figure_3

%% Figure 3 
%  rank histogram

subplot(2,1,1) 
fname     = 'obs_diag_output.nc';
timeindex = -1;
plotdat   = plot_rank_histogram(fname, timeindex, 'AMSRE_BRIGHTNESS_T');


