% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% Select field to plot (U, V, W, GZ, T, MU, QV, QC, QR)

field_num = input('Input field type, 1=U, 2=V, 3=W, 4=GZ, 5=T, 6=MU, 7=QV, 8=QC, 9=QR, 10=XLAND: ');

map_proj = {'lambert', 'ups', 'mercator'};
state_name = {'True_State' 'Prior_Diag' 'Posterior_Diag'};

disp('Plotting horizontal map of a linear combinaison of'); disp (state_name)

fact = input('Input coefficients (between [ ]):');

if fact(1) ~= 0.0
   fname = char(state_name(1));
elseif fact(2) ~= 0.0
   fname = char(state_name(2));
elseif fact(3) ~= 0.0
   fname = char(state_name(3));
else
   error('Nothing to plot.')
end

xlon = getnc(fname, 'XLON');
we = size(xlon, 2);
xlat = getnc(fname, 'XLAT');
sn = size(xlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);

cop = zeros(1,3);
cop(1) = -1;

if (fact(2) ~= 0.0) | (fact(3) ~= 0.0)
   if fact(2) ~= 0.0
      fname = char(state_name(2));
   else
      fname = char(state_name(3));
   end
   ncopy = getnc(fname, 'copy');
   ens_size = size(ncopy, 1) - 2;
   disp(['The first ' int2str(ens_size) ' copies are ensemble members.'])
   disp(['The ensemble mean is copy #' int2str(ens_size+1)])
   disp(['The ensemble spread is copy #',int2str(ens_size+2)])
   cop(2:3) = input('Input copy for Prior and/or Posterior: ');
end

%--Read data
nc = netcdf( [fname, '.nc'] , 'read' ) ;
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
if (field_num == 6) | (field_num == 10)
   field_level = 1;
else
   field_level = input('Input level: ');
end

nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
var_name = 'U';
iso = [0.5:1:5];
if field_num > 1
   nx = we;
   ny = sn + 1;
   var_units = 'V (m/s)';
   var_name = 'V';
end
if field_num > 2
   nx = we;
   ny = sn;
   var_units = 'W (m/s)';
   var_name = 'W';
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
   iso = [0.5:0.5:5];
end
if field_num > 5
   var_units = 'MU (Pa)';
   var_name = 'MU';
   iso = [100:100:600];
end
if field_num > 6
   var_units = 'QV (kg/kg)';
   var_name = 'QVAPOR';
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
if field_num > 9
   var_units = 'XLAND (-)';
   var_name = 'XLAND';
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.8*scrsz(4) 0.8*scrsz(4)])

m = ceil(sqrt(ftime-stime+1));

pane = 1;

for itime = stime:ftime

plot_title = [var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)];

% Extract field

field = zeros(sn,we);

for istate = 1:3

if fact(istate) ~= 0.0

if (field_num ~= 6) & (field_num < 10)
   corner = [itime cop(istate) field_level -1 -1];
   end_point = [itime cop(istate) field_level -1 -1];
   stride = [1 1 1 1 1];
end
if field_num == 6
   corner = [itime cop(istate) -1 -1];
   end_point = [itime cop(istate) -1 -1];
   stride = [1 1 1 1];
end
if field_num == 10
   corner = [-1 -1];
   end_point = [-1 -1];
   stride = [1 1];
end

   fname = char(state_name(istate));
   stag_field = getnc(fname, var_name,corner,end_point,stride);
   if field_num == 1
      for iy = 1:sn
      for ix = 1:we
         field(iy,ix) = field(iy,ix) + fact(istate)*(stag_field(iy,ix) + stag_field(iy,ix+1))/2.0;
      end
      end
   elseif field_num == 2
      for iy = 1:sn
      for ix = 1:we
         field(iy,ix) = field(iy,ix) + fact(istate)*(stag_field(iy,ix) + stag_field(iy+1,ix))/2.0;
      end
      end
   else
      field = field + fact(istate)*stag_field;
   end

end

end

% Plot field

subplot(m,m,pane);

axesm(map_proj{mp},'Origin',[0 cen_lon 0],'MapParallels',[stdlat1 stdlat2],...
      'MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]); framem;

plotm(coast,'color',[0 0 0]);
plotm(usalo('statebvec'),'color',[0 0 0]);
plotm(usalo('conusvec'),'color',[0 0 0]);

% axis( [-0.6 0.6 .2 1.2 ]) % This works pretty well for present CONUS domain

if min(min(field)) ~= max(max(field))

[C h] = contourm(xlat,xlon,field);
h = clabelm(C,h,'labelspacing',288);  set(h,'Fontsize',9);

title(plot_title)

%[C,h] = contour (field, iso);
%hold on
%[Cm,hm] = contour (field, -iso, '--');
%colorbar('vert')
%clabel(C, h);
%clabel(Cm, hm);

end

pane = pane + 1;

end

% Loop for another try
%map_wrf_diff_time;
