% psfc_map

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

field_name = 'PSFC';

map_proj = {'lambert', 'ups', 'mercator'};

fname = 'psfc';

nc = netcdf( [fname,'.nc'] , 'nowrite' ) ;

stdlat1 = nc.TRUELAT1(:);
stdlat2 = nc.TRUELAT2(:);
cen_lon = nc.CEN_LON(:);
mp = nc.MAP_PROJ(:);
dt = nc.DT(:);

close(nc);

lon = getnc(fname, 'XLONG');
xlon = squeeze(lon(1,:,:));
we = size(xlon, 2);
lat = getnc(fname, 'XLAT');
xlat = squeeze(lat(1,:,:));
sn = size(xlat, 1);

f_size = we*sn;

for iy = 1:sn
   for ix = 1:we
      if(xlon(iy,ix) > 0.0)
         xlon(iy,ix) = xlon(iy,ix) - 360.0;
      end
   end
end

minlat = min(xlat(:)); maxlat = max(xlat(:));
minlon = min(xlon(:)); maxlon = max(xlon(:));

true_times = getnc(fname, 'Times');
num_true_times = size(true_times, 1)

stime = input('Initial time : ');
ftime = input('End time : ');

var_units = ' (Pa)';
iso = [-15:2:15];

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.8*scrsz(4) 0.8*scrsz(4)])

m = ceil(sqrt(ftime-stime+1));

pane = 1;

x = [1:ftime-stime+1];
rmse = x;

for itime = stime:ftime

   plot_title = [ field_name var_units ...
			    '   ' true_times(itime,:) ];

% Extract field

   field1 = getnc(fname, field_name,[itime -1 -1],[itime -1 -1],[1 1 1]);
   field2 = getnc(fname, field_name,[itime+1 -1 -1],[itime+1 -1 -1],[1 1 1]);

field = (field2 - field1)/dt;

field_vec = reshape(field,f_size,1);

rmse(pane) = sqrt((field_vec'*field_vec)/(f_size));

% Plot field

   subplot(m,m,pane);

   axesm(map_proj{mp},'Origin',[0 cen_lon 0],'MapParallels', ...
	 [stdlat1 stdlat2],...
	 'MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]);
   framem;

   plotm(coast,'color',[0 0 0]);
   plotm(usalo('statebvec'),'color',[0 0 0]);
   plotm(usalo('conusvec'),'color',[0 0 0]);

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
