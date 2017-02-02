%% psfc_movie

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

field_name = 'MU';

fname = 'psfc.nc';

times    = nc_varget(fname,     'time'); Ntimes = size(times, 1);
xlon     = nc_varget(fname, 'XLON_d01');  we    = size( xlon, 2);
xlat     = nc_varget(fname, 'XLAT_d01');  sn    = size( xlat, 1);
level    = nc_varget(fname,'level_d01');  bt    = size(level, 1);
dt       = nc_varget(fname,       'DT');

stime = input('Initial time : ');
ftime = input('End time : ');

var_units = ' (Pa)';
iso       = -2:0.1:2;

%set(gca,'nextplot','replacechildren');

pane = 1;

for itime = stime:ftime

   %% Extract two adjacent timesteps of the field and manipulate

   start = [itime    1  1] -1;
   count = [    2   -1 -1];

   field12 = nc_varget(fname, field_name, start, count);
   field   = (field12(2,:,:) - field12(1,:,:))/dt;

   % Plot field

   if min(min(field)) ~= max(max(field))

      [C h] = contourf(field, iso); caxis([min(iso(:)),max(iso(:))]);
      %  h = clabel(C,h,'labelspacing',288);  set(h,'Fontsize',12);
      %  hold on
      %  [Cm hm] = contourm(xlat,xlon,field, -iso, 'b--','LineWidth',2);
      %  hm = clabelm(Cm,hm,'labelspacing',288);  set(hm,'Fontsize',12);

      plot_title = [ field_name var_units '   ' times(itime,:) ];

      title(plot_title)
      colorbar

      F(pane) = getframe(gcf);

      pane = pane + 1;

   end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
