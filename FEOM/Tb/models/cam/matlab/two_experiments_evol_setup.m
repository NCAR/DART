% This is the one

SE_itwin = 'SE_RUD_itwin1_obs_diag.nc';
FV_itwin = 'FV_RUD_itwin2_obs_diag.nc';
% FV_ftwin = 'FV_RUD_ftwin1_obs_diag.nc';
% SE_ftwin = 'SE_RUD_ftwin1_obs_diag.nc';
% FV_adap1 = 'FV_RUD_itwin_ainf1_obs_diag.nc';
% SE_adap1 = 'SE_RUD_ftwin_ainf1_obs_diag.nc';
% FV_adap2 = 'FV_RUD_ftwin_ainf2_obs_diag.nc';

% plotdat = plot_rmse_xxx_evolution(FV_adap2, 'bias');
% plotdat = plot_rmse_xxx_evolution(FV_adap2, 'spread');
% plotdat = plot_rmse_xxx_profile(  FV_adap2, 'bias');
% plotdat = plot_rmse_xxx_profile(  FV_adap2, 'spread');

% files    = {FV_itwin, SE_itwin, FV_ftwin, SE_ftwin, FV_adap2};
files    = { FV_itwin,   SE_itwin};
titles   = {'FV identical twin', 'SE identical twin'};
varnames = {'RADIOSONDE_V_WIND_COMPONENT', 'RADIOSONDE_TEMPERATURE'};
% varnames = {'RADIOSONDE_TEMPERATURE'};
qtty     = 'rmse';
prpo     = 'forecast';
% Not CAM levels; obs_diag levels.
%    plevel = 1000.,850.,700.,500.,400.,300.,200.,150.,100.,40.,15.,5.

levelind = 1;

two_experiments_evolution(files, titles, varnames, qtty, prpo, levelind)

% Plot the evolution of the spread.

% copystring = 'spread';
% plotdat    = plot_evolution(FV_ftwin, 'totalspread','RADIOSONDE_TEMPERATURE');
% 
% foreach FILE ( RADIOSONDE_* )
%   mv $FILE FV_$FILE
% end
% 
% varnames = { 'RADIOSONDE_U_WIND_COMPONENT'};
% varnames = { 'RADIOSONDE_V_WIND_COMPONENT'}; 
% 
% print -f1 -dpdf SE_FV_rmse_NH_500hPa.pdf
% print -f2 -dpdf SE_FV_rmse_SH_500hPa.pdf
% print -f3 -dpdf SE_FV_rmse_TR_500hPa.pdf
% print -f4 -dpdf SE_FV_rmse_NA_500hPa.pdf


