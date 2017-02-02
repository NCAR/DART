%% psfc_map

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

field_name = 'PSFC';

map_proj = {'lambert', 'ups', 'mercator'};

fname = 'psfc.nc';

%% Read map metadata
times    = nc_varget(fname,     'time'); Ntimes = size(times, 1);
xlon     = nc_varget(fname, 'XLON_d01');  we    = size( xlon, 2);
xlat     = nc_varget(fname, 'XLAT_d01');  sn    = size( xlat, 1);
level    = nc_varget(fname,'level_d01');  bt    = size(level, 1);
stdlat1  = nc_varget(fname, 'TRUELAT1');
stdlat2  = nc_varget(fname, 'TRUELAT2');
cen_lat  = nc_varget(fname,  'CEN_LAT');
cen_lon  = nc_varget(fname,  'CEN_LON');
dt       = nc_varget(fname,       'DT');
mp       = nc_varget(fname, 'MAP_PROJ')+1;  % convert [0,N-1] to [1,N]
map_proj = nc_attget(fname, 'MAP_PROJ','units');

true_times     = nc_varget(fname, 'Times');
num_true_times = size(true_times, 1);

stime = input('Initial time : ');
ftime = input('End time : ');

f_size = we*sn;

minlat = min(xlat(:)); maxlat = max(xlat(:));
minlon = min(xlon(:)); maxlon = max(xlon(:));

var_units = ' (Pa)';
iso       = -15:2:15;

scrsz     = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.8*scrsz(4) 0.8*scrsz(4)])

m = ceil(sqrt(ftime-stime+1));

pane = 1;

x = 1:ftime-stime+1;
rmse = x;

for itime = stime:ftime

   plot_title = [ field_name var_units '   ' true_times(itime,:) ];

   % Extract two adjacent timesteps of the field and manipulate

   start   = [itime    1  1] -1;
   count   = [    2   -1 -1];
   field12 = nc_varget(fname, field_name, start, count);
   field   = (field12(2,:,:) - field12(1,:,:))/dt;

   field_vec = reshape(field,f_size,1);

   rmse(pane) = sqrt((field_vec'*field_vec)/(f_size));

   % Plot field

   subplot(m,m,pane);

   axesm(map_proj{mp},'Origin',[0 cen_lon 0],'MapParallels', ...
	 [stdlat1 stdlat2], ...
	 'MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]);
   framem;

   plotm(coast,             'color',[0 0 0]);
   plotm(usalo('statebvec'),'color',[0 0 0]);
   plotm(usalo( 'conusvec'),'color',[0 0 0]);

   %axis( [-0.65 0.65 .1 1.45 ]) % This works pretty well for present CONUS domain

   if min(min(field)) ~= max(max(field))

      %     [C h] = contourm(xlat,xlon,field, iso, 'r','LineWidth',2);
      [C h] = contourfm(xlat,xlon,field, iso); caxis([min(iso(:)),max(iso(:))]);
      %     h = clabelm(C,h,'labelspacing',288);  set(h,'Fontsize',12);
      hold on
      %     [Cm hm] = contourm(xlat,xlon,field, -iso, 'b--','LineWidth',2);
      %     hm = clabelm(Cm,hm,'labelspacing',288);  set(hm,'Fontsize',12);

   end

   title(plot_title)
   colorbar

   pane = pane + 1;

end

figure(2);
plot(x,rmse,'LineWidth',2)

% Loop for another try
%map_wrf_diff_time;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
