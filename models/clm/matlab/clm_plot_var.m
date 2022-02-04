function clm_plot_var(x)
%% DART clm_plot_var - plots the data structure that results from clm_get_var.
%
% EXAMPLE 1:
% fname      = 'clm_preassim_mean.2000-01-05-00000.nc';
% varname    = 'frac_sno';
% levelindex = 1;
% timeindex  = 1;
% x = clm_get_var(fname,varname,levelindex,timeindex);
% clm_plot_var(x);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%

h = imagesc(x.lonarray,x.latarray,x.datmat);
set(h,'AlphaData',~isnan(x.datmat))
set(gca, 'YDir', 'normal', 'FontSize', 18)
continents('hollow')
colorbar
str = sprintf('%s level %d time index %d',x.varname,x.levelindex,x.timeindex);
ht = title(str);
set(ht,'Interpreter','none');
hx = xlabel(x.filename);
set(hx,'Interpreter','none','FontSize',12);

