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
fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_Mg_MOPnIAS_20M_100km_p10p30/obs_diag_output.nc';
%
copystring1    = 'ens_mean';
copystring2    = 'observation';
copystring3    = 'bias';
copystring4    = 'rmse';
%
%plot1 = plot_profile(fname,copystring1);
%plot2 = plot_profile(fname,copystring2);
plot3 = plot_profile(fname,copystring3);
plot4 = plot_profile(fname,copystring4);
