function obsstruct = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
                                     QCString, maxQC, verbose)
%
% fname         = 'obs_sequence_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'NCEP BUFR observation';
% QCString      = 'DART quality control';
% maxQC         = 2;
% verbose       = 1;   % anything > 0 == 'true'
%
% bob = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, maxQC, verbose);
%
% view(0,90);   % for a traditional '2D' plot

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

% Read the observation sequence

obsstruct = read_obs_netcdf(fname, ObsTypeString, region, ...
                    CopyString, QCString, maxQC, verbose);

%-------------------------------------------------------------------------------
% Create graphic with area-weighted symbols for the good observations.
%-------------------------------------------------------------------------------

figure(1); clf

xmin = min(region(1:2));
xmax = max(region(1:2));
ymin = min(region(3:4));
ymax = max(region(3:4));
zmin = min(obsstruct.z);
zmax = max(obsstruct.z);

scalearray = scaleme(obsstruct.obs, 36);
scatter3(obsstruct.lons, obsstruct.lats, obsstruct.z, ...
         scalearray, obsstruct.obs, 'd', 'filled');

axis([xmin xmax ymin ymax zmin zmax])

str1 = sprintf('%s level (%.2f - %.2f)',ObsTypeString,zmin,zmax);
str2 = sprintf('%s (%d locations)',CopyString,length(obsstruct.obs));
str3 = sprintf('%s - %s',obsstruct.timestring(1,:),obsstruct.timestring(2,:));

title( {str1, str3, str2}, 'Interpreter','none','FontSize',16);
xlabel('longitude')
ylabel('latitude')

if     (obsstruct.Ztyp(1) == -2) % VERTISUNDEF     = -2
   zlabel('curious ... undefined')
elseif (obsstruct.Ztyp(1) == -1) % VERTISSURFACE   = -1
   zlabel('surface')
elseif (obsstruct.Ztyp(1) ==  1) % VERTISLEVEL     =  1
   zlabel('level')
elseif (obsstruct.Ztyp(1) ==  2) % VERTISPRESSURE  =  2
   set(gca,'ZDir','reverse')
   zlabel('pressure')
elseif (obsstruct.Ztyp(1) ==  3) % VERTISHEIGHT    =  3
   zlabel('height')
end

myworldmap;
set(gca,'CLim',[min(obsstruct.obs) max(obsstruct.obs)])
h = colorbar;
set(get(h,'YLabel'),'String',ObsTypeString,'Interpreter','none')

%-------------------------------------------------------------------------------
% Create graphic of spatial distribution of 'bad' observations & their QC value.
%-------------------------------------------------------------------------------

if (obsstruct.numbadqc > 0 ) 

   figure(2); clf
   
   subplot('position',[0.1 0.20 0.8 0.65])
   scalearray = 128 * ones(size(obsstruct.badobs.obs));
   
   zmin = min(obsstruct.badobs.z);
   zmax = max(obsstruct.badobs.z);
   
   scatter3(obsstruct.badobs.lons, obsstruct.badobs.lats, obsstruct.badobs.z, ...
            scalearray, obsstruct.badobs.qc,'filled')

   axis([xmin xmax ymin ymax zmin zmax])
   
   str2 = sprintf('%s (%d bad observations)',CopyString,length(obsstruct.badobs.obs));
   
   title( {str1, str3, str2}, 'Interpreter','none','FontSize',16);
   xlabel('longitude')
   ylabel('latitude')
   
   if     (obsstruct.badobs.Ztyp(1) == -2) % VERTISUNDEF     = -2
      zlabel('curious ... undefined')
   elseif (obsstruct.badobs.Ztyp(1) == -1) % VERTISSURFACE   = -1
      zlabel('surface')
   elseif (obsstruct.badobs.Ztyp(1) ==  1) % VERTISLEVEL     =  1
      zlabel('level')
   elseif (obsstruct.badobs.Ztyp(1) ==  2) % VERTISPRESSURE  =  2
      set(gca,'ZDir','reverse')
      zlabel('pressure')
   elseif (obsstruct.badobs.Ztyp(1) ==  3) % VERTISHEIGHT    =  3
      zlabel('height')
   end
   
   myworldmap;
   set(gca,'CLim',[min(obsstruct.badobs.qc) max(obsstruct.badobs.qc)])
   h = colorbar;
   set(get(h,'YLabel'),'String',QCString,'Interpreter','none')
   
   subplot('position',[0.1 0.05 0.8 0.10])
   axis off
   
   qcvals  = unique(obsstruct.badobs.qc);
   qccount = zeros(size(qcvals));
   for i = 1:length(qcvals)
      qccount(i) = sum(obsstruct.badobs.qc == qcvals(i));
      s{i} = sprintf('%d obs with qc == %d',qccount(i),qcvals(i));
   end
   
   dy = 1.0/length(s);
   for i = 1:length(s)
      text(0.0, (i-1)*dy ,s{i})
   end

end





function h = myworldmap

%---------------------------------------------------------------------------
% GET THE ELEVATION DATA AND SET UP THE ASSOCIATED COORDINATE DATA
%---------------------------------------------------------------------------

load topo;                 % GET Matlab-native [180x360] ELEVATION DATASET
lats = [-89.5:89.5];       % CREATE LAT ARRAY FOR TOPO MATRIX
lons = [0.5:359.5];        % CREATE LON ARRAY FOR TOPO MATRIX
nlon = length(lons);
nlat = length(lats);

%---------------------------------------------------------------------------
% IF WE NEED TO SWAP HEMISPHERES, DO SO NOW.
% If we didn't explicitly tell it, make a guess.
%---------------------------------------------------------------------------

ax   = axis;

if (ax(1) < -2)
   lons = lons - 180.0;
   topo = [ topo(:,nlon/2+1:nlon) topo(:,1:nlon/2) ];
end

%---------------------------------------------------------------------------
% We need to determine the geographic subset of the elevation matrix.
%---------------------------------------------------------------------------

lon_ind1 = min(find(ax(1) <= lons));
lon_ind2 = min(find(ax(2) <= lons));
lat_ind1 = min(find(ax(3) <= lats));
lat_ind2 = min(find(ax(4) <= lats));

if (isempty(lon_ind1)) lon_ind1 = 1;    end;
if (isempty(lon_ind2)) lon_ind2 = nlon; end;
if (isempty(lat_ind1)) lat_ind1 = 1;    end;
if (isempty(lat_ind2)) lat_ind2 = nlat; end;

elev = topo(lat_ind1:lat_ind2,lon_ind1:lon_ind2);
x    = lons(lon_ind1:lon_ind2);
y    = lats(lat_ind1:lat_ind2);

%---------------------------------------------------------------------------
% Contour the "subset"
% There are differences between 6.5 and 7.0 that make changing the colors
% of the filled contours a real pain. Providing both solutions.
%---------------------------------------------------------------------------

orgholdstate = ishold;
hold on;

switch  get(gca,'ZDir')
   case 'reverse'
      zlevel = max(ax(5:6));
   otherwise
      zlevel = min(ax(5:6));
end

fcolor = [0.7 0.7 0.7];    % light grey

[c,h] = contourf(x,y,elev,[0.0 0.0],'k-');

new_level = 1000;
h_patch   = get(h, 'Children');

for i = 1:numel(h_patch)
    y = get(h_patch(i), 'YData');
    s = size(y);
    set(h_patch(i), 'ZData', zlevel*ones(s),'FaceColor',fcolor);
end

if (orgholdstate == 0) hold off; end;




function s = scaleme(x,minsize)
% scaleme returns a uniformly scaled array the same size as the input
% array where the maximum is 10 times the minimum 
maxsize = 10*minsize;
minx    = min(x);
maxx    = max(x);
slope   = (maxsize-minsize)/(maxx-minx);
b       = minsize - slope*minx;

s = x*slope + b;

