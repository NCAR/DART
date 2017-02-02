function clm_plot_var(x)
%% DART clm_plot_var - plots the data structure that results from clm_get_var.
%
% EXAMPLE 1:
% fname      = 'Prior_Diag.2000-01-05-00000.nc';
% varname    = 'frac_sno';
% levelindex = 1;
% timeindex  = 1;
% copystring = 'ensemble member 3';
% x = clm_get_var(fname,varname,copystring,levelindex,timeindex);
% clm_plot_var(x);
%
% EXAMPLE 2: as above, compare to comparable field in CLM history file.
% clmfname   = '/glade/scratch/thoar/enstest_0907/enstest_0907.clm2_0003.r.2000-01-05-00000.nc';
% x = clm_get_var(fname,varname,copystring,levelindex,timeindex,clmfname);
% clm_plot_var(x);
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (isfield(x,'clmres'))
   subplot(2,1,1)
end
   h1 = imagesc(x.lon, x.lat, x.datmat);
   set(h1,'AlphaData',~isnan(x.datmat))
   set(gca,'YDir','normal')
   h2 = title(sprintf('DART version of %s',x.varname));
   set(h2,'Interpreter','none')
   worldmap;
   colorbar;

if (isfield(x,'clmres'))
   subplot(2,1,2)
      h3 = imagesc(x.lon, x.lat, x.clmres);
      set(h3,'AlphaData',~isnan(x.clmres))
      set(gca,'YDir','normal')
      h4 = title(sprintf('Restart file version of %s',x.varname));
      set(h4,'Interpreter','none')
      worldmap;
      colorbar;
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
