fname         = '/glade/scratch/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_60M_NL/obs_diag_output.nc';
%
timeindex     = 19;
plot_var      = 'RADIOSONDE_TEMPERATURE';
plot_var      = 'MOPITT_CO_RETRIEVAL';
plot          = plot_rank_histogram(fname, timeindex, plot_var);
