fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_1_Trc_Ret/obs_epoch_001.nc';
region        = [0 360 -90 90 -Inf Inf];
ObsTypeString = 'MOPITT_CO_RETRIEVAL';
CopyString1    = 'NCEP BUFR observation';
CopyString2    = 'posterior ensemble mean';
QCString      = 'DART quality control';
maxgoodQC     = 2;
verbose       = 1;
twoup         = 0;
plot          = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
                      QCString, maxgoodQC, verbose, twoup);
