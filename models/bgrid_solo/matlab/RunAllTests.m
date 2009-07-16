function RunAllTests(dummy)
% RunAllTests.m

%------------------------------------------------------------
% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$
%------------------------------------------------------------

if (nargin() > 0)
   interactive = 1;
else
   interactive = 0;
end

if (interactive)
 plot_bins
 fprintf('Finished %s pausing, hit any key\n','plot_bins'); pause
 plot_ens_err_spread
 fprintf('Finished %s pausing, hit any key\n','plot_ens_err_spread'); pause
 plot_ens_time_series
 fprintf('Finished %s pausing, hit any key\n','plot_ens_time_series'); pause
 plot_ens_mean_time_series
 fprintf('Finished %s pausing, hit any key\n','plot_ens_mean_time_series'); pause
end

 clear pinfo; close all; 
 pinfo.truth_file  = 'True_State.nc';
 pinfo.diagn_file  = 'Prior_Diag.nc';
 pinfo.truth_time  = [1 2];
 pinfo.diagn_time  = [1 2];
 pinfo.model       = 'FMS_Bgrid';
 pinfo.fname       = 'Prior_Diag.nc';
 pinfo.var         = 'u';
 pinfo.level       = 2;
 pinfo.levelindex  = 2;
 pinfo.longitude   = 258.00;
 pinfo.lonindex    = 43;
 pinfo.latitude    = 42.00;
 pinfo.latindex    = 22;

 close all; PlotBins(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotBins'); pause

 clf; PlotEnsErrSpread(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsErrSpread'); pause

 clf; PlotEnsTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsTimeSeries'); pause

 clf; PlotEnsMeanTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsMeanTimeSeries'); pause

%------------------------------------------------------------
%plot_correl
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_correl
 fprintf('Finished %s pausing, hit any key\n','plot_correl'); pause
end

 clear pinfo; clf
 pinfo.model              = 'FMS_Bgrid';
 pinfo.num_state_vars     = 4;
 pinfo.num_ens_members    = 24;
 pinfo.time_series_length = 2;
 pinfo.min_ens_mem        = 1;
 pinfo.max_ens_mem        = 24;
 pinfo.vars               = {'ps'  't'  'u'  'v'};
 pinfo.fname              = 'Prior_Diag.nc';
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
 fprintf('Finished %s pausing, hit any key\n','PlotCorrel'); pause

%------------------------------------------------------------
%plot_phase_space
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_phase_space
 fprintf('Finished %s pausing, hit any key\n','plot_phase_space'); pause
end

 clear pinfo; clf
 pinfo.fname       = 'True_State.nc';
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
 pinfo.ens_mem     = 'true state';
 pinfo.ltype       = 'k-';
 
 PlotPhaseSpace(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotPhaseSpace'); pause

%------------------------------------------------------------
%plot_reg_factor
%------------------------------------------------------------
% plot_reg_factor

%------------------------------------------------------------
%plot_sawtooth
%------------------------------------------------------------
if (interactive)
 clear; close all; plot_sawtooth
 fprintf('Finished %s pausing, hit any key\n','plot_sawtooth'); pause
end

 clear pinfo; close all
 pinfo.model              = 'FMS_Bgrid';
 pinfo.vars               = {'ps'  't'  'u'  'v'};
 pinfo.truth_file         = 'True_State.nc';
 pinfo.prior_file         = 'Prior_Diag.nc';
 pinfo.posterior_file     = 'Posterior_Diag.nc';
 pinfo.var_names          = 'u';
 pinfo.level              = 3;
 pinfo.levelindex         = 3;
 pinfo.latitude           = 42.0000;
 pinfo.longitude          = 258.0000;
 pinfo.latindex           = 22;
 pinfo.lonindex           = 43;
 pinfo.copies             = 3;
 pinfo.copyindices        = [7 12 17];
 pinfo.truth_time         = [1 2];
 pinfo.prior_times        = [1 2];
 pinfo.posterior_times    = [1 2];

 PlotSawtooth(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotSawtooth'); pause

%------------------------------------------------------------
%plot_smoother_err
%------------------------------------------------------------
% plot_smoother_err

%------------------------------------------------------------
%plot_total_err
%------------------------------------------------------------
if (interactive)
 clear; close all; plot_total_err
 fprintf('Finished %s pausing, hit any key\n','plot_total_err'); pause
end

 clear pinfo; close all
 pinfo.model              = 'FMS_Bgrid';
 pinfo.num_state_vars     = 4;
 pinfo.num_ens_members    = 24;
 pinfo.time_series_length = 2;
 pinfo.min_ens_mem        = 1;
 pinfo.max_ens_mem        = 24;
 pinfo.vars               = {'ps'  't'  'u'  'v'};
 pinfo.truth_file         = 'True_State.nc';
 pinfo.diagn_file         = 'Prior_Diag.nc';
 pinfo.truth_time         = [1 2];
 pinfo.diagn_time         = [1 2];

 PlotTotalErr(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotTotalErr'); pause

%------------------------------------------------------------
%plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_var_var_correl
 fprintf('Finished %s pausing, hit any key\n','plot_var_var_correl'); pause
end

 clear pinfo; clf
 pinfo.fname       = 'Prior_Diag.nc';
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
 fprintf('Finished %s pausing, hit any key\n','PlotVarVarCorrel'); pause

%------------------------------------------------------------
%plot_jeff_correl - virtually identical to plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_jeff_correl
 fprintf('Finished %s pausing, hit any key\n','plot_jeff_correl'); pause
end

 PlotJeffCorrel(pinfo)

