fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_60M_p20p20/obs_epoch_001.nc';
region        = [0 360 -90 90 -Inf Inf];
ObsTypeString = 'MOPITT_CO_RETRIEVAL';
CopyString    = 'NCEP BUFR observation';
QCString      = 'DART quality control';
maxgoodQC     = 4;
verbose       = 1;
twoup         = 2;
plot          = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
                      QCString, maxgoodQC, verbose, twoup);
