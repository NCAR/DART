function RunAllTests(dummy)
%% RunAllTests.m bgrid_solo

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

 pinfo = CheckModelCompatibility('true_state.nc','preassim.nc');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.diagn_file);
 pinfo.var         = 'u';
 pinfo.level       = 2;
 pinfo.levelindex  = 2;
 pinfo.longitude   = 258.00;
 pinfo.lonindex    = 43;
 pinfo.latitude    = 42.00;
 pinfo.latindex    = 22;
 pinfo.fname       = pinfo.diagn_file;

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

 pinfo = CheckModel('preassim.nc');
 pinfo.base_var           = 't';
 pinfo.comp_var           = 'u';
 pinfo.base_time          = 0;
 pinfo.base_tmeind        = 1;
 pinfo.base_lvl           = 2;
 pinfo.base_lvlind        = 2;
 pinfo.base_lat           = 39.0000;
 pinfo.base_latind        = 22;
 pinfo.base_lon           = 255.0000;
 pinfo.base_lonind        = 43;
 pinfo.comp_lvl           = 2;
 pinfo.comp_lvlind        = 2;

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

 pinfo.fname       = 'true_state.nc';
 pinfo.model       = 'FMS_Bgrid';
 pinfo.var1name    = 'ps';
 pinfo.var2name    = 't';
 pinfo.var3name    = 'u';
 
 pinfo.var1_lvl    = 1;
 pinfo.var2_lvl    = 1;
 pinfo.var3_lvl    = 1;
 pinfo.var1_lvlind = 1;
 pinfo.var2_lvlind = 1;
 pinfo.var3_lvlind = 1;
 pinfo.var1_lat    = 39.0000;
 pinfo.var2_lat    = 39.0000;
 pinfo.var3_lat    = 36.0000;
 pinfo.var1_latind = 22;
 pinfo.var2_latind = 22;
 pinfo.var3_latind = 21;
 pinfo.var1_lon    = 255.0000;
 pinfo.var2_lon    = 255.0000;
 pinfo.var3_lon    = 252.0000;
 pinfo.var1_lonind = 43;
 pinfo.var2_lonind = 43;
 pinfo.var3_lonind = 42;
 pinfo.ens_mem     = 'ensemble member 1';
 pinfo.ltype       = 'k-';
 
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

 pinfo    = CheckModelCompatibility('preassim.nc','analysis.nc');
 pinfo.prior_time     = pinfo.truth_time;
 pinfo.prior_file     = pinfo.truth_file;
 pinfo.posterior_time = pinfo.diagn_time;
 pinfo.posterior_file = pinfo.diagn_file;
 pinfo.truth_file     = 'true_state.nc';
 pinfo = rmfield(pinfo,{'diagn_file','truth_time','diagn_time'});
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.prior_file);
 pinfo.var_names      = 'ps';
 pinfo.level          = 1;
 pinfo.levelindex     = 1;
 pinfo.latitude       = 39.0;
 pinfo.latindex       = 22;
 pinfo.longitude      = 255.0;
 pinfo.lonindex       = 43;
 pinfo.copyindices    = [5 10 15];
 pinfo.copies         = length(pinfo.copyindices);

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

 pinfo    = CheckModelCompatibility('true_state.nc','preassim.nc');

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

 pinfo  = CheckModel('preassim.nc');
 pinfo.model       = 'FMS_Bgrid';
 pinfo.base_var    = 't';
 pinfo.comp_var    = 'u';
 pinfo.base_time   = 0;
 pinfo.base_tmeind = 1;
 pinfo.base_lvl    = 4;
 pinfo.base_lat    = -39.0000;
 pinfo.base_lon    = 183.0000;
 pinfo.base_latind = 9;
 pinfo.base_lvlind = 4;
 pinfo.base_lonind = 31;
 pinfo.comp_lvl    = 2;
 pinfo.comp_lat    = 42.0000;
 pinfo.comp_lon    = 258.0000;
 pinfo.comp_lvlind = 2;
 pinfo.comp_latind = 22;
 pinfo.comp_lonind = 43;

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
