function RunAllTests(dummy)
%% RunAllTests.m pe2lyr

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
 pinfo.var        = 'u';
 pinfo.levelindex = 1;
 pinfo.lonindex   = 69;
 pinfo.latindex   = 14;
 pinfo.level      = 1;
 pinfo.longitude  = 255.0;
 pinfo.latitude   = 38.966;
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
 pinfo.base_var           = 'v';
 pinfo.comp_var           = 'u';
 pinfo.base_tmeind        =   2;
 pinfo.base_lvlind        =   1;
 pinfo.base_lonind        =  69;
 pinfo.base_latind        =  14;
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvl           =     1;
 pinfo.base_lon           = 255.0;
 pinfo.base_lat           =  38.96;
 pinfo.comp_lvlind        =    2;
 pinfo.comp_lvl           =    2;
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
 pinfo.var1name     = 'u';
 pinfo.var2name     = 'v';
 pinfo.var3name     = 'z';
 pinfo.var1_lvlind  = 1;
 pinfo.var2_lvlind  = 1;
 pinfo.var3_lvlind  = 1;
 pinfo.var1_lvl     = 1;
 pinfo.var2_lvl     = 1;
 pinfo.var3_lvl     = 1;
 pinfo.var1_latind  = 14;
 pinfo.var2_latind  = 14;
 pinfo.var3_latind  = 14;
 pinfo.var1_lat     = 38.96;
 pinfo.var2_lat     = 38.96;
 pinfo.var3_lat     = 38.96;
 pinfo.var1_lonind  = 69;
 pinfo.var2_lonind  = 69;
 pinfo.var3_lonind  = 69;
 pinfo.var1_lon     = 255;
 pinfo.var2_lon     = 255;
 pinfo.var3_lon     = 255;
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
 pinfo.var_names  = 'z';
 pinfo.levelindex  = 1;
 pinfo.latindex    = 14;
 pinfo.lonindex    = 69;
 pinfo.level       = 1;
 pinfo.latitude    = 38.96;
 pinfo.longitude   = 255;
 pinfo.copies      = 3;
 pinfo.copyindices = [5 10 15];

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
 vars1 = rmfield(vars1,{'time','time_series_length','fname'});
 vars2 = CheckModelCompatibility(truth_file,diagn_file);
 pinfo = CombineStructs(vars1,vars2);
 pinfo.levels = ncread(pinfo.diagn_file,'lev');
 pinfo.lons   = ncread(pinfo.diagn_file,'lon');
 pinfo.lats   = ncread(pinfo.diagn_file,'lat');

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
 pinfo.base_var           = 'u';
 pinfo.comp_var           = 'v';
 pinfo.base_tmeind        = 2;
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvl           = 1;
 pinfo.base_lvlind        = 1;
 pinfo.base_lat           = 38.966;
 pinfo.base_latind        = 14;
 pinfo.base_lon           = 255.0;
 pinfo.base_lonind        = 69;
 pinfo.comp_lvl           = 1;
 pinfo.comp_lvlind        = 1;
 pinfo.comp_lat           = -38.966;
 pinfo.comp_latind        = 35;
 pinfo.comp_lon           = pinfo.base_lon;
 pinfo.comp_lonind        = pinfo.base_lonind;

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
 PlotJeffCorrel(pinfo)
 fprintf('Finished %s\n','PlotJeffCorrel')

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
