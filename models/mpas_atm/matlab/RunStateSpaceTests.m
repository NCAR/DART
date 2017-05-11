function RunStateSpaceTests(dummy)
%% RunStateSpaceTests.m

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

if (interactive)
 fprintf('Starting %s\n','plot_bins');
 plot_bins
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

 clear pinfo; close all;

 truth_file = 'true_state.nc';
 diagn_file = 'preassim.nc';
 vars1 = CheckModel(diagn_file);
 vars1 = rmfield(vars1,{'time','time_series_length','fname'});
 vars2 = CheckModelCompatibility(truth_file,diagn_file);
 pinfo = CombineStructs(vars1,vars2);
 pinfo.var        = 'theta';
 pinfo.lonCell    = nc_varget(truth_file,'lonCell');
 pinfo.latCell    = nc_varget(truth_file,'latCell');
 pinfo.area       = nc_varget(truth_file,'areaCell');
 pinfo.lonunits   = 'degrees_east';
 pinfo.latunits   = 'degrees_north';
 pinfo.levelindex = 1;
 pinfo.cellindex  = 5121;
 pinfo.level      = 21;
 pinfo.longitude  = pinfo.lonCell(pinfo.cellindex);
 pinfo.latitude   = pinfo.latCell(pinfo.cellindex);
 clear vars1 vars2

 fprintf('Starting %s\n','PlotBins');
 close all; PlotBins(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotBins'); pause

 fprintf('Starting %s\n','PlotEnsErrSpread');
 close all; PlotEnsErrSpread(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsErrSpread'); pause

 fprintf('Starting %s\n','PlotEnsTimeSeries');
 close all; PlotEnsTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsTimeSeries'); pause

 fprintf('Starting %s\n','PlotEnsMeanTimeSeries');
 close all; PlotEnsMeanTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsMeanTimeSeries'); pause

%------------------------------------------------------------
%plot_correl
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_correl');
 clear; close all; plot_correl
 fprintf('Finished %s pausing, hit any key\n','plot_correl'); pause
end

 clear pinfo;
 pinfo                    = CheckModel('preassim.nc');
 pinfo.lonCell            = nc_varget(pinfo.fname,'lonCell');
 pinfo.latCell            = nc_varget(pinfo.fname,'latCell');
 pinfo.area               = nc_varget(pinfo.fname,'areaCell');
 pinfo.lonunits           = 'degrees_east';
 pinfo.latunits           = 'degrees_north';
 pinfo.base_var           = 'uReconstructZonal';
 pinfo.comp_var           = 'uReconstructZonal';
 pinfo.base_tmeind        =  1;
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvlind        = 21;
 pinfo.base_lvl           = 21;
 pinfo.base_cellindex     = 5121;
 pinfo.base_lon           = pinfo.lonCell(pinfo.base_cellindex);
 pinfo.base_lat           = pinfo.latCell(pinfo.base_cellindex);
 pinfo.comp_lvlind        = 21;
 pinfo.comp_lvl           = 21;
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);

 fprintf('Starting %s\n','PlotCorrel');
 PlotCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotCorrel'); pause

%------------------------------------------------------------
%plot_phase_space
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_phase_space');
 clear; close all; plot_phase_space
 fprintf('Finished %s pausing, hit any key\n','plot_phase_space'); pause
end

 clear pinfo;
 pinfo                 = CheckModel('preassim.nc');
 pinfo.lonCell         = nc_varget(pinfo.fname,'lonCell');
 pinfo.latCell         = nc_varget(pinfo.fname,'latCell');
 pinfo.area            = nc_varget(pinfo.fname,'areaCell');
 pinfo.lonunits        = 'degrees_east';
 pinfo.latunits        = 'degrees_north';
 pinfo.var1name        = 'theta';
 pinfo.var2name        = 'w';
 pinfo.var3name        = 'rho';
 pinfo.var1_lvl        = 21;
 pinfo.var1_lvlind     = 21;
 pinfo.var1_cellindex  = 5121;
 pinfo.var1_lat        = pinfo.latCell(pinfo.var1_cellindex);
 pinfo.var1_lon        = pinfo.lonCell(pinfo.var1_cellindex);
 pinfo.var2_lvl        = 21;
 pinfo.var2_lvlind     = 21;
 pinfo.var2_cellindex  = 5121;
 pinfo.var2_lat        = pinfo.latCell(pinfo.var2_cellindex);
 pinfo.var2_lon        = pinfo.lonCell(pinfo.var2_cellindex);
 pinfo.var3_lvl        = 21;
 pinfo.var3_lvlind     = 21;
 pinfo.var3_cellindex  = 5121;
 pinfo.var3_lat        = pinfo.latCell(pinfo.var3_cellindex);
 pinfo.var3_lon        = pinfo.lonCell(pinfo.var3_cellindex);
 pinfo.ens_mem         = 'ensemble mean';
 pinfo.ltype           = 'k-';

 fprintf('Starting %s\n','PlotPhaseSpace');
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
 fprintf('Starting %s\n','plot_sawtooth');
 clear; close all; plot_sawtooth
 fprintf('Finished %s pausing, hit any key\n','plot_sawtooth'); pause
end

 clear pinfo; clf
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
 pinfo.lonCell        = nc_varget(truth_file,'lonCell');
 pinfo.latCell        = nc_varget(truth_file,'latCell');
 pinfo.var_names  = 'theta';
 pinfo.levelindex  = 21;
 pinfo.level       = 21;
 pinfo.cellindex   = 5121;
 pinfo.latitude    = pinfo.latCell(pinfo.cellindex);
 pinfo.longitude   = pinfo.lonCell(pinfo.cellindex);
 pinfo.copies      = 1;
 pinfo.copyindices = [ 1 ];

 fprintf('Starting %s\n','PlotSawtooth');
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
 fprintf('Starting %s\n','plot_total_err');
 clear; close all; plot_total_err
 fprintf('Finished %s pausing, hit any key\n','plot_total_err'); pause
end

 clear pinfo; clf
 truth_file           = 'true_state.nc';
 diagn_file           = 'preassim.nc';
 vars1                = CheckModel(diagn_file);
 rmfield(vars1,{'time','time_series_length','fname'});
 vars2                = CheckModelCompatibility(truth_file,diagn_file);
 pinfo                = CombineStructs(vars1,vars2);
 pinfo.num_state_vars = 2; % just do the first two
 pinfo.area           = nc_varget(truth_file,'areaCell');

 fprintf('Starting %s\n','PlotTotalErr');
 PlotTotalErr(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotTotalErr'); pause

%------------------------------------------------------------
%plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_var_var_correl');
 clear; close all; plot_var_var_correl
 fprintf('Finished %s pausing, hit any key\n','plot_var_var_correl'); pause
end

% Darwin, Australia = -12.435235, 130.841675
% Tahiti            : -17.643368, 210.575638

 clear pinfo; clf
 diagn_file = 'preassim.nc';
 pinfo = CheckModel(diagn_file);
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 pinfo.lonCell            = nc_varget(pinfo.fname,'lonCell');
 pinfo.latCell            = nc_varget(pinfo.fname,'latCell');
 pinfo.base_var           = 'surface_pressure';
 pinfo.comp_var           = 'surface_pressure';
 pinfo.base_tmeind        = 2;
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvl           = 1;
 pinfo.base_lvlind        = 1;
 pinfo.base_cellindex     = 5680;
 pinfo.base_lat           = pinfo.latCell(pinfo.base_cellindex);
 pinfo.base_lon           = pinfo.lonCell(pinfo.base_cellindex);
 pinfo.comp_lvl           = 1;
 pinfo.comp_lvlind        = 1;
 pinfo.comp_cellindex     = 1193;
 pinfo.comp_lat           = pinfo.latCell(pinfo.comp_cellindex);
 pinfo.comp_lon           = pinfo.lonCell(pinfo.comp_cellindex);

 fprintf('Starting %s\n','PlotVarVarCorrel');
 PlotVarVarCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotVarVarCorrel'); pause

%------------------------------------------------------------
%plot_jeff_correl - correlation evolution
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_jeff_correl');
 clear; close all; plot_jeff_correl
 fprintf('Finished %s pausing, hit any key\n','plot_jeff_correl'); pause
end

 figure(1); clf

 % Largely uses same pinfo as PlotVarVarCorrel
 pinfo = rmfield(pinfo,'base_tmeind');

 fprintf('Starting %s\n','PlotJeffCorrel');
 PlotJeffCorrel(pinfo)
 fprintf('Finished %s\n','PlotJeffCorrel')

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
