function RunAllTests(dummy)
%% RunAllTests.m wrf

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin() > 0)
   interactive = 1;
else
   interactive = 0;
end

%%

figure(1)
if (interactive)
 fprintf('Starting %s\n','plot_bins');
 clear truth_file diagn_file; close all; plot_bins
 fprintf('Finished %s ... pausing, hit any key\n','plot_bins'); pause
 fprintf('Starting %s\n','plot_ens_err_spread');
 plot_ens_err_spread
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_err_spread'); pause
 fprintf('Starting %s\n','plot_ens_time_series');
 plot_ens_time_series
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_time_series'); pause
 fprintf('Starting %s\n','plot_ens_mean_time_series');
 plot_ens_mean_time_series
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_mean_time_series'); pause
end

 fprintf('Starting %s\n','PlotBins');
 clear pinfo; close all;

 truth_file = 'true_state.nc';
 diagn_file = 'preassim.nc';
 vars1 = CheckModel(diagn_file);
 vars1 = rmfield(vars1,{'time','time_series_length','fname'});
 vars2 = CheckModelCompatibility(truth_file,diagn_file);
 pinfo = CombineStructs(vars1,vars2);
 pinfo.lonmat     = ncread(pinfo.truth_file, 'XLONG_d01');
 pinfo.latmat     = ncread(pinfo.truth_file, 'XLAT_d01');
 [nlats, nlons]   = size(pinfo.lonmat);
 pinfo.var        = 'T_d01';
 pinfo.levelindex = 1;
 pinfo.lonindex   = nlons;
 pinfo.latindex   = nlats;
 pinfo.level      = 2;
 pinfo.longitude  = pinfo.lonmat(pinfo.latindex,pinfo.lonindex);
 pinfo.latitude   = pinfo.latmat(pinfo.latindex,pinfo.lonindex);
 pinfo.fname      = pinfo.diagn_file;
 clear vars1 vars2
 
 PlotBins(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotBins'); pause

 fprintf('Starting %s\n','PlotEnsErrSpread');
 close all; PlotEnsErrSpread(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotEnsErrSpread'); pause

 fprintf('Starting %s\n','PlotEnsTimeSeries');
 close all; PlotEnsTimeSeries(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotEnsTimeSeries'); pause

 fprintf('Starting %s\n','PlotEnsMeanTimeSeries');
 close all; PlotEnsMeanTimeSeries(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotEnsMeanTimeSeries'); pause

%% ----------------------------------------------------------
% plot_correl
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n', 'plot_correl');
 clear diagn_file; close all; plot_correl
 fprintf('Finished %s ... pausing, hit any key\n','plot_correl'); pause
end

 fprintf('Starting %s\n','PlotCorrel');
 clear pinfo; clf

 pinfo                    = CheckModel('preassim.nc');
 pinfo.time               = ncread(pinfo.fname,'time');
 pinfo.time_series_length = length(pinfo.time);
 pinfo.lonmat             = ncread(pinfo.fname, 'XLONG_d01');
 pinfo.latmat             = ncread(pinfo.fname, 'XLAT_d01');
 [nlats, nlons]           = size(pinfo.lonmat);
 pinfo.base_var           = 'T_d01';
 pinfo.comp_var           = 'T_d01';
 pinfo.base_tmeind        =   3;
 pinfo.base_lvlind        =  40;
 pinfo.base_lonind        =  round(nlons/2);
 pinfo.base_latind        =  round(nlats/2);
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvl           = pinfo.base_lvlind;
 pinfo.base_lon           = pinfo.lonmat(pinfo.base_latind,pinfo.base_lonind);
 pinfo.base_lat           = pinfo.latmat(pinfo.base_latind,pinfo.base_lonind);
 pinfo.comp_lvlind        = pinfo.base_lvlind;
 pinfo.comp_lvl           = pinfo.base_lvl;
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);

 PlotCorrel(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotCorrel'); pause

%% -----------------------------------------------------------
% plot_phase_space
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_phase_space');
 clear fname; close all; plot_phase_space
 fprintf('Finished %s ... pausing, hit any key\n','plot_phase_space'); pause
end

 fprintf('Starting %s\n','PlotPhaseSpace');
 clear pinfo; clf

 pinfo              = CheckModel('preassim.nc');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 pinfo.var1name     = 'T_d01';
 pinfo.var2name     = 'PH_d01';
 pinfo.var3name     = 'W_d01';
 pinfo.lonmat       = ncread(pinfo.fname, 'XLONG_d01');
 pinfo.latmat       = ncread(pinfo.fname, 'XLAT_d01');
 [nlats, nlons]     = size(pinfo.lonmat);
 pinfo.var1_lvlind  = 1;
 pinfo.var2_lvlind  = 1;
 pinfo.var3_lvlind  = 1;
 pinfo.var1_lvl     = 1;
 pinfo.var2_lvl     = 1;
 pinfo.var3_lvl     = 1;
 pinfo.var1_latind  = round(nlats/2);
 pinfo.var2_latind  = round(nlats/2);
 pinfo.var3_latind  = round(nlats/2);
 pinfo.var1_lonind  = round(nlons/2);
 pinfo.var2_lonind  = round(nlons/2);
 pinfo.var3_lonind  = round(nlons/2);
 pinfo.var1_lat     = pinfo.latmat(pinfo.var1_latind, pinfo.var1_lonind);
 pinfo.var2_lat     = pinfo.latmat(pinfo.var2_latind, pinfo.var2_lonind);
 pinfo.var3_lat     = pinfo.latmat(pinfo.var3_latind, pinfo.var3_lonind);
 pinfo.var1_lon     = pinfo.lonmat(pinfo.var1_latind, pinfo.var1_lonind);
 pinfo.var2_lon     = pinfo.lonmat(pinfo.var2_latind, pinfo.var2_lonind);
 pinfo.var3_lon     = pinfo.lonmat(pinfo.var3_latind, pinfo.var3_lonind);
 pinfo.ens_mem      = 'ensemble mean';
 pinfo.ltype        = 'k-';

 PlotPhaseSpace(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotPhaseSpace'); pause

%% ----------------------------------------------------------
% plot_reg_factor
%------------------------------------------------------------
% plot_reg_factor

%% ----------------------------------------------------------
% plot_sawtooth
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_sawtooth');
 clear truth_file posterior_file prior_file; close all; plot_sawtooth
 fprintf('Finished %s ... pausing, hit any key\n','plot_sawtooth'); pause
end

 fprintf('Starting %s\n','PlotSawtooth');
 clear pinfo; close all

 truth_file     = 'true_state.nc';
 prior_file     = 'preassim.nc';
 posterior_file = 'analysis.nc';
 pinfo = CheckModelCompatibility(prior_file,posterior_file);
 pinfo.prior_time     = pinfo.truth_time;
 pinfo.posterior_time = pinfo.diagn_time;
 pinfo.truth_file     = truth_file;
 pinfo.prior_file     = prior_file;
 pinfo.posterior_file = posterior_file;
 pinfo = rmfield(pinfo,{'diagn_file','truth_time','diagn_time'});
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(prior_file);
 pinfo.var_names      = 'U_d01';
 pinfo.latmat         = ncread(pinfo.prior_file, 'XLAT_U_d01');
 pinfo.lonmat         = ncread(pinfo.prior_file, 'XLONG_U_d01');
[nlats, nlons]        = size(pinfo.lonmat);
 pinfo.levelindex     = 1;
 pinfo.latindex       = round(nlats/2);
 pinfo.lonindex       = round(nlons/2);
 pinfo.level          = 1;
 pinfo.latitude       = pinfo.latmat(pinfo.latindex, pinfo.lonindex);
 pinfo.longitude      = pinfo.lonmat(pinfo.latindex, pinfo.lonindex);
 pinfo.copies         = 3;
 pinfo.copyindices    = [1 10 20];

 PlotSawtooth(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotSawtooth'); pause



%% ----------------------------------------------------------
% plot_total_err
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_total_err');
 clear truth_file diagn_file; close all; plot_total_err
 fprintf('Finished %s ... pausing, hit any key\n','plot_total_err'); pause
end

 fprintf('Starting %s\n','PlotTotalErr');
 clear pinfo; clf

 truth_file = 'true_state.nc';
 diagn_file = 'preassim.nc';
 vars1 = CheckModel(diagn_file);
 rmfield(vars1,{'time','time_series_length','fname'});
 vars2 = CheckModelCompatibility(truth_file,diagn_file);
 pinfo = CombineStructs(vars1,vars2);
 pinfo.num_state_vars = 2; % just do the first two

 PlotTotalErr(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotTotalErr'); pause

%% ----------------------------------------------------------
% plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_var_var_correl');
 clear fname; close all; plot_var_var_correl
 fprintf('Finished %s ... pausing, hit any key\n','plot_var_var_correl'); pause
end

 fprintf('Starting %s\n','PlotVarVarCorrel');
 clear pinfo; clf

 diagn_file = 'preassim.nc';
 pinfo = CheckModel(diagn_file);
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 pinfo.base_var     = 'T_d01';
 pinfo.comp_var     = 'PH_d01';
 pinfo.base_tmeind  = 3;
 pinfo.base_time    = pinfo.time(pinfo.base_tmeind);
 pinfo.lonmat       = ncread(pinfo.fname, 'XLONG_d01');
 pinfo.latmat       = ncread(pinfo.fname, 'XLAT_d01');
 [nlats, nlons]     = size(pinfo.lonmat);
 pinfo.base_lvlind  = 1;
 pinfo.base_lvl     = 1;
 pinfo.base_latind  = round(nlats/3);
 pinfo.base_lonind  = round(nlons/2);
 pinfo.base_lat     = pinfo.latmat(pinfo.base_latind, pinfo.base_lonind);
 pinfo.base_lon     = pinfo.lonmat(pinfo.base_latind, pinfo.base_lonind);
 pinfo.comp_lvlind  = 2;
 pinfo.comp_lvl     = 2;
 pinfo.comp_latind  = round(nlats/6);
 pinfo.comp_lonind  = round(nlons/2);
 pinfo.comp_lat     = pinfo.latmat(pinfo.comp_latind, pinfo.comp_lonind);
 pinfo.comp_lon     = pinfo.lonmat(pinfo.comp_latind, pinfo.comp_lonind);

 PlotVarVarCorrel(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotVarVarCorrel'); pause

%% ----------------------------------------------------------
% plot_jeff_correl - correlation evolution
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_jeff_correl');
 clear fname; close all; plot_jeff_correl
 fprintf('Finished %s ... pausing, hit any key\n','plot_jeff_correl'); pause
end

 fprintf('Starting %s\n','PlotJeffCorrel');
 clf
 PlotJeffCorrel(pinfo)
 fprintf('Finished %s\n','PlotJeffCorrel')

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
