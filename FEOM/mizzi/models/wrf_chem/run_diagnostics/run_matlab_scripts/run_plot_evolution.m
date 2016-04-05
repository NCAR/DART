%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_40M_p30p30_sp4/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_40M_p30p30_sp4/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_MgDA_40M_p30p30_sp4/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p10p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p20p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p30p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p30p30/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_MgDA_20M_p10p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p10p30/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_bar_1_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_bar_2_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_loc_a_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_loc_b_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_bar_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_loc_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_MgDA_20M_100km_p10p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_loc_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_bloc_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_bloc_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p75/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p50/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p60/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_DBL_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_DBL_NV_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_DBL_bloc_p10p30/obs_diag_output.nc';
fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_Mg_MOPnIAS_20M_100km_p10p30/obs_diag_output.nc';
%
npar=1;
copystring    = {'spread'};
copystring    = {'totalspread'};
copystring    = {'rmse'};
%
nvar=2;
varname      = {'MOPITT_CO_RETRIEVAL','IASI_CO_RETRIEVAL'};
for ipar=1:npar
for ivar=1:nvar
    plot = plot_evolution(fname,copystring{ipar},varname{ivar});
end
end
