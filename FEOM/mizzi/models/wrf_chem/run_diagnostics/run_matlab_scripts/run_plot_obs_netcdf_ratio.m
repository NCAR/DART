%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_TRet_20M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_TRet_20M_HL/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_RtDA_20M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_RtDA_20M_HL/obs_epoch_019.nc';
%
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_TRet_60M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_TRet_60M_HL/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_RtDA_60M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_RtDA_60M_HL/obs_epoch_019.nc';
%
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_TRet_20M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_TRet_20M_HL/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_20M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_20M_HL/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_TRet_60M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_TRet_60M_HL/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_60M/obs_epoch_019.nc';
fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_60M_HL/obs_epoch_019.nc';
%
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_TRet_20M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_TRet_20M_HL/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_RtDA_20M/obs_epoch_019.nc';
%fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_RtDA_20M_HL/obs_epoch_019.nc';
%
region        = [0 360 -90 90 -Inf Inf];
ObsTypeString = 'MOPITT_CO_RETRIEVAL';
CopyString1    = 'prior ensemble spread';
CopyString2    = 'observation error variance';
QCString      = 'DART quality control';
maxgoodQC     = 2;
verbose       = 1;
twoup         = 0;
plot          = plot_obs_netcdf_ratio(fname, ObsTypeString, region, CopyString1, ...
				      CopyString2, QCString, maxgoodQC, verbose, twoup);
