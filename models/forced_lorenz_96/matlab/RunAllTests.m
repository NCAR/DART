function RunAllTests(dummy)
%% RunAllTests.m

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

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
 pinfo.truth_file     = './True_State.nc';
 pinfo.diagn_file     = './Prior_Diag.nc';
 pinfo.model          = 'Forced_Lorenz_96';
 pinfo.var            = 'state';
 pinfo.truth_time     = [1 1000];
 pinfo.diagn_time     = [1 1000];
 pinfo.var_inds       = [1 13 27];
 
 clf; PlotBins(pinfo)
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
 pinfo.model              =  'Forced_Lorenz_96';
 pinfo.def_var            =  'state';
 pinfo.num_state_vars     =  80;
 pinfo.num_model_vars     =  40;
 pinfo.num_force_vars     =  40;
 pinfo.num_ens_members    =  24;
 pinfo.time_series_length =  1000;
 pinfo.min_state_var      =  1;
 pinfo.max_state_var      =  80;
 pinfo.min_model_var      =  1;
 pinfo.max_model_var      =  40;
 pinfo.min_force_var      =  1;
 pinfo.max_force_var      =  40;
 pinfo.min_ens_mem        =  1;
 pinfo.max_ens_mem        =  24;
 pinfo.def_state_vars     =  [1 13 27];
 pinfo.def_force_vars     =  [41 53 67];
 pinfo.fname              =  'Prior_Diag.nc';
 pinfo.base_var           =  'state';
 pinfo.base_var_index     =  30;
 pinfo.base_time          =  300;
 
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
 pinfo.fname    = 'True_State.nc';
 pinfo.model    = 'Forced_Lorenz_96';
 pinfo.var1name = 'state';
 pinfo.var2name = 'state';
 pinfo.var3name = 'state';
 pinfo.var1ind  = 10;
 pinfo.var2ind  = 20;
 pinfo.var3ind  = 30;
 pinfo.ens_mem  = 'true state';
 pinfo.ltype    = 'k-';

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
 pinfo.model              = 'Forced_Lorenz_96';
 pinfo.def_var            = 'state';
 pinfo.def_state_vars     = [1 13 27];
 pinfo.def_force_vars     = [41 53 67];
 pinfo.prior_file         = 'Prior_Diag.nc';
 pinfo.posterior_file     = 'Posterior_Diag.nc';
 pinfo.diagn_file         = 'Prior_Diag.nc';
 pinfo.diagn_time         = [1 -1];
 pinfo.truth_file         = 'True_State.nc';
 pinfo.truth_time         = [1 -1];
 pinfo.var                = 'state';
 pinfo.var_inds           = [1 20 40];
 pinfo.copyindices        = [7 12 17];

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
 pinfo.model              = 'Forced_Lorenz_96';
 pinfo.def_var            = 'state';
 pinfo.def_state_vars     = [1 13 27];
 pinfo.def_force_vars     = [41 53 67];
 pinfo.truth_file         = 'True_State.nc';
 pinfo.diagn_file         = 'Prior_Diag.nc';
 pinfo.truth_time         = [1 1000];
 pinfo.diagn_time         = [1 1000];

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
 pinfo.fname           = 'Prior_Diag.nc';
 pinfo.model           = 'Forced_Lorenz_96';
 pinfo.base_var        = 'state';
 pinfo.state_var       = 'state';
 pinfo.base_var_index  = 30;
 pinfo.base_time       = 300;
 pinfo.state_var_index = 60;

 PlotVarVarCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotVarVarCorrel'); pause

%------------------------------------------------------------
%plot_jeff_correl - virtually identical to plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_jeff_correl
 fprintf('Finished %s pausing, hit any key\n','plot_jeff_correl'); pause
end

 clear pinfo; clf
 pinfo.fname           = 'Prior_Diag.nc';
 pinfo.base_var        = 'state';
 pinfo.state_var       = 'state';
 pinfo.base_var_index  = 20;
 pinfo.base_time       = 300;
 pinfo.state_var_index = 10;

 PlotJeffCorrel(pinfo)

