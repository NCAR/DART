% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Select field to plot (U, V, W, GZ, T, MU, QV, QC, QR)

field_num = input('Input field type, 1=U, 2=V, 3=W, 4=GZ, 5=T, 6=MU, 7=QV, 8=QC, 9=QR: ');

map_proj = {'lambert', 'ups', 'mercator'};

% Get file name of true state file
fname = 'True_State';
xlon = getnc(fname, 'XLON');
we = size(xlon, 2);
xlat = getnc(fname, 'XLAT');
sn = size(xlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);

%--Read data
nc = netcdf( 'True_State.nc' , 'read' ) ;
stdlat1 = nc.TRUELAT1(:);
stdlat2 = nc.TRUELAT2(:);
cen_lat = nc.CEN_LAT(:);
cen_lon = nc.CEN_LON(:);
mp = nc.MAP_PROJ(:);

close(nc)

minlat = min(xlat(:)); maxlat = max(xlat(:));
minlon = min(xlon(:)); maxlon = max(xlon(:));

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1)

stime = input('Initial time : ');
ftime = input('End time : ');

% Get level for free atmosphere fields
if field_num == 6
   field_level = 1;
else
   field_level = input('Input level: ');
end

var_units = 'U (m/s)';
var_name = 'U';
maxlev = bt;
iso = [0.5:1:5];
if field_num > 1
var_units = 'V (m/s)';
   var_name = 'V';
end
if field_num > 2
var_units = 'W (m/s)';
   var_name = 'W';
   maxlev = bt + 1;
   iso = [0.01:0.01:0.1];
end
if field_num > 3
var_units = 'GZ (m^2/s^2)';
   var_name = 'PH';
   iso = [50:50:300];
end
if field_num > 4
var_units = 'T (K)';
   var_name = 'T';
   maxlev = bt;
   iso = [0.5:0.5:5];
end
if field_num > 5
var_units = 'MU (Pa)';
   var_name = 'MU';
   maxlev = 1;
   iso = [100:100:600];
end
if field_num > 6
var_units = 'QV (kg/kg)';
   var_name = 'QVAPOR';
   maxlev = bt;
   iso = [0.0001:0.0001:0.001];
end
if field_num > 7
var_units = 'QC (kg/kg)';
   var_name = 'QCLOUD';
   iso = [0.00001:0.00001:0.0001];
end
if field_num > 8
var_units = 'QR (kg/kg)';
   var_name = 'QRAIN';
   iso = [0.00001:0.00001:0.0001];
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.8*scrsz(4) 0.8*scrsz(4)])

m = ceil(sqrt(ftime-stime+1));

pane = 1;

for itime = stime:ftime

plot_title = ['True state   ' var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)];

if maxlev > 1
   stride = [1 1 1 1 1];
   corner_t = [itime -1 field_level -1 -1];
   end_point_t = [itime -1 field_level -1 -1];
else
   stride = [1 1 1 1];
   corner_t = [itime -1 -1 -1];
   end_point_t = [itime -1 -1 -1];
end

% Extract field

stag_field = getnc(fname, var_name,corner_t,end_point_t,stride);

field = zeros(sn,we);
if field_num == 1
for iy = 1:sn
for ix = 1:we
   field(iy,ix) = (stag_field(iy,ix) + stag_field(iy,ix+1))/2.0;
end
end
elseif field_num == 2
for iy = 1:sn
for ix = 1:we
   field(iy,ix) = (stag_field(iy,ix) + stag_field(iy+1,ix))/2.0;
end
end
else
field = stag_field;
end

% Plot field

subplot(m,m,pane);

axesm(map_proj{mp},'Origin',[0 cen_lon 0],'MapParallels',[stdlat1 stdlat2],...
      'MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]); framem;

plotm(coast,'color',[0 0 0]);
plotm(usalo('statebvec'),'color',[0 0 0]);
plotm(usalo('conusvec'),'color',[0 0 0]);

% axis( [-0.6 0.6 .2 1.2 ]) % This works pretty well for present CONUS domain

if min(field) ~= max(field)

[C h]=contourm(xlat,xlon,field) ;
h = clabelm(C,h,'labelspacing',288);  set(h,'Fontsize',9);

% Can get lat,lon lines via, e.g., contourm(xlat,xlon,xlat) 
[C h]=contourm(xlat,xlon,xlat,'k') ;
h = clabelm(C,h,'labelspacing',288);  set(h,'Fontsize',9);

title(plot_title)
%colorbar('vert')

end

pane = pane + 1;

end

% Loop for another try
%map_wrf;
