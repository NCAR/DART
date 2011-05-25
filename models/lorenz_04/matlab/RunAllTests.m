function RunAllTests(dummy)
%% RunAllTests.m

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
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
 pinfo.model          = 'Lorenz_04';
 pinfo.var            = 'state';
 pinfo.truth_time     = [1 100];
 pinfo.diagn_time     = [1 100];
 pinfo.var_inds       = [100 200 300];
 
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
pinfo.model              = 'Lorenz_04';
pinfo.def_var            = 'state';
pinfo.num_state_vars     = 960;
pinfo.num_ens_members    = 24;
pinfo.time_series_length = 100;
pinfo.min_state_var      = 1;
pinfo.max_state_var      = 960;
pinfo.min_ens_mem        = 1;
pinfo.max_ens_mem        = 24;
pinfo.def_state_vars     = [1 320 640];
pinfo.fname              = 'Prior_Diag.nc';
pinfo.base_var           = 'state';
pinfo.base_var_index     = 100;
pinfo.base_time          = 20;
 
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
 pinfo.model    = 'Lorenz_04';
 pinfo.var1name = 'state';
 pinfo.var2name = 'state';
 pinfo.var3name = 'state';
 pinfo.var1ind  = 100;
 pinfo.var2ind  = 200;
 pinfo.var3ind  = 300;
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
 pinfo.model           = 'Lorenz_04';
 pinfo.def_var         = 'state';
 pinfo.def_state_vars  = [1 320 640];
 pinfo.prior_file      = 'Prior_Diag.nc';
 pinfo.posterior_file  = 'Posterior_Diag.nc';
 pinfo.truth_file      = 'True_State.nc';
 pinfo.truth_time      = [1 -1];
 pinfo.var             = 'state';
 pinfo.var_inds        = [1 320 640];
 pinfo.copyindices     = [7 12 17];
 pinfo.prior_times     = [1 100];
 pinfo.posterior_times = [1 100];

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
 pinfo.def_var            = 'state';
 pinfo.num_state_vars     = 960;
 pinfo.num_ens_members    = 24;
 pinfo.time_series_length = 100;
 pinfo.min_state_var      = 1;
 pinfo.max_state_var      = 960;
 pinfo.min_ens_mem        = 1;
 pinfo.max_ens_mem        = 24;
 pinfo.def_state_vars     = [1 320 640];
 pinfo.truth_file         = 'True_State.nc';
 pinfo.diagn_file         = 'Prior_Diag.nc';
 pinfo.truth_time         = [1 100];
 pinfo.diagn_time         = [1 100];

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
 pinfo.model           = 'Lorenz_04';
 pinfo.base_var        = 'state';
 pinfo.state_var       = 'state';
 pinfo.base_var_index  = 200;
 pinfo.base_time       = 30;
 pinfo.state_var_index = 800;

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
 pinfo.base_var_index  = 800;
 pinfo.base_time       = 30;
 pinfo.state_var_index = 200;

 PlotJeffCorrel(pinfo)

