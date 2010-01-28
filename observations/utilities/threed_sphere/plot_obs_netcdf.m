function obsstruct = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
                                     QCString, maxQC, verbose)
%% plot_obs_netcdf will plot the locations and values of the observations in a DART netcdf file.
%     any observations with a QC value greater than 'maxgoodQC' will get
%     plotted on a separate figure ... color-coded to its QC value, not the
%     observation value.
%
%--------------------------------------------------
% EXAMPLE 1: plotting just one type of observation
%--------------------------------------------------
% fname         = 'obs_sequence_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'NCEP BUFR observation';
% QCString      = 'DART quality control';
% maxgoodQC     = 2;
% verbose       = 1;   % anything > 0 == 'true'
%
% bob = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, maxgoodQC, verbose);
%
% view(0,90);   % for a traditional '2D' plot
%
%--------------------------------------------------
% EXAMPLE 2: plotting all the observation types 
%--------------------------------------------------
% fname         = 'obs_sequence_001.nc';
% ObsTypeString = 'ALL';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'WOD observation';
% QCString      = 'WOD QC';
% maxgoodQC     = 0;
% verbose       = 1;   % anything > 0 == 'true'
%
% bob = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, maxgoodQC, verbose);

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
%  <next few lines under version control, do not edit>
%  $URL$
%  $Id$
%  $Revision$
%  $Date$

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

%% Read the observation sequence

obsstruct = read_obs_netcdf(fname, ObsTypeString, region, ...
                    CopyString, QCString, maxQC, verbose);

%% Create graphic with area-weighted symbols for the good observations.
%  It has happened that there have been zero good observations in a file.

xmin = min(region(1:2));
xmax = max(region(1:2));
ymin = min(region(3:4));
ymax = max(region(3:4));
zmin = min(obsstruct.z);
zmax = max(obsstruct.z);

pstruct.colorbarstring = obsstruct.ObsTypeString;
pstruct.region = region;
pstruct.str3   = sprintf('%s - %s',obsstruct.timestring(1,:),obsstruct.timestring(2,:));

if ( length(obsstruct.obs) < 1 ) 
   fprintf('There are no ''good'' observations to plot\n')
else

   figure(gcf+1); clf

   % choose a symbol size based on the number of obs to plot.

   if (length(obsstruct.obs) < 1000) 
      pstruct.scalearray = scaleme(obsstruct.obs, 36);
   else
      pstruct.scalearray = 50.0 * ones(size(obsstruct.obs));
   end
   pstruct.clim   = [min(obsstruct.obs) max(obsstruct.obs)];
   pstruct.str2   = sprintf('%s (%d locations)',obsstruct.CopyString,length(obsstruct.obs));

   % If all the observations live on the same level ... make a 2D plot.

   if ( zmin ~= zmax )

      pstruct.axis = [xmin xmax ymin ymax zmin zmax];
      pstruct.str1 = sprintf('%s level (%.2f - %.2f)',obsstruct.ObsTypeString,zmin,zmax);

      plot_3D(obsstruct, pstruct)

   else

      pstruct.axis = [xmin xmax ymin ymax];
      pstruct.str1 = sprintf('%s',obsstruct.ObsTypeString);

      plot_2D(obsstruct, pstruct)

   end
end

%% Create graphic of spatial distribution of 'bad' observations & their QC value.
%
% 0     observation assimilated
% 1     observation evaluated only
%   --- everything above this means the prior and posterior are OK
% 2     assimilated, but the posterior forward operator failed
% 3     Evaluated only, but the posterior forward operator failed
%   --- everything above this means only the prior is OK
% 4     prior forward operator failed
% 5     not used
% 6     prior QC rejected
% 7     outlier rejected

dartqc_strings = { ...
   '''observation evaluated only''', ...
   '''assimilated, but the posterior forward operator failed''', ...
   '''evaluated only, but the posterior forward operator failed''',...
   '''prior forward operator failed''',...
   '''not used''',...
   '''prior QC rejected''',...
   '''outlier rejected''',...
   '''reserved for future use'''};

if (obsstruct.numbadqc > 0 ) % if there are bad observation to plot ... carry on.

   figure(gcf+1); clf
   
   subplot('position',[0.1 0.20 0.8 0.65])
   
   zmin = min(obsstruct.badobs.z);
   zmax = max(obsstruct.badobs.z);

   pstruct.scalearray = 128 * ones(size(obsstruct.badobs.obs));
   pstruct.colorbarstring = QCString;
   pstruct.clim = [min(obsstruct.badobs.qc) max(obsstruct.badobs.qc)];
   pstruct.str1 = sprintf('%s level (%.2f - %.2f)',obsstruct.ObsTypeString,zmin,zmax);
   pstruct.str2 = sprintf('%s (%d bad observations)',  ...
                                    obsstruct.CopyString, ...
                             length(obsstruct.badobs.obs));
 
   obsstruct.badobs.obs = obsstruct.badobs.qc;  % plot QC values, not obs values
   if ( zmin ~= zmax )

      pstruct.axis = [xmin xmax ymin ymax zmin zmax];

      plot_3D(obsstruct.badobs, pstruct)

   else

      pstruct.axis = [xmin xmax ymin ymax];

      plot_2D(obsstruct.badobs, pstruct)

   end
   
   subplot('position',[0.1 0.05 0.8 0.10])
   axis off

   %% If the QC is a DART QC, we know how to interpret them.

   switch lower(strtrim(QCString))
      case 'dart quality control',

         qcvals  = unique(obsstruct.badobs.qc); 
         qccount = zeros(size(qcvals));
         for i = 1:length(qcvals)
            qccount(i) = sum(obsstruct.badobs.qc == qcvals(i));
            s{i} = sprintf('%d obs with qc == %d %s',qccount(i),qcvals(i), ...
                   dartqc_strings{qcvals(i)});
         end
   
         dy =  1.0/length(s);
         for i = 1:length(s)
            text(0.0, (i-1)*dy ,s{i})
         end

      otherwise,
         str = sprintf('no way to interpret values of %s',strtrim(QCString));
         text(0.0, 0.0, str)
   end
end



function h = myworldmap

%%--------------------------------------------------------------------------
% GET THE ELEVATION DATA AND SET UP THE ASSOCIATED COORDINATE DATA
%---------------------------------------------------------------------------

load topo;               % GET Matlab-native [180x360] ELEVATION DATASET
lats = -89.5:89.5;       % CREATE LAT ARRAY FOR TOPO MATRIX
lons = 0.5:359.5;        % CREATE LON ARRAY FOR TOPO MATRIX
nlon = length(lons);
nlat = length(lats);

%%--------------------------------------------------------------------------
% IF WE NEED TO SWAP HEMISPHERES, DO SO NOW.
% If we didn't explicitly tell it, make a guess.
%---------------------------------------------------------------------------

ax   = axis;

if (ax(1) < -2)
   lons = lons - 180.0;
   topo = [ topo(:,nlon/2+1:nlon) topo(:,1:nlon/2) ];
end

%%--------------------------------------------------------------------------
% We need to determine the geographic subset of the elevation matrix.
%---------------------------------------------------------------------------

lon_ind1 = min(find(ax(1) <= lons));
lon_ind2 = min(find(ax(2) <= lons));
lat_ind1 = min(find(ax(3) <= lats));
lat_ind2 = min(find(ax(4) <= lats));

if (isempty(lon_ind1)), lon_ind1 = 1;    end;
if (isempty(lon_ind2)), lon_ind2 = nlon; end;
if (isempty(lat_ind1)), lat_ind1 = 1;    end;
if (isempty(lat_ind2)), lat_ind2 = nlat; end;

elev = topo(lat_ind1:lat_ind2,lon_ind1:lon_ind2);
x    = lons(lon_ind1:lon_ind2);
y    = lats(lat_ind1:lat_ind2);

%%--------------------------------------------------------------------------
% Contour the "subset"
% There are differences between 6.5 and 7.0 that make changing the colors
% of the filled contours a real pain.
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

h_patch   = get(h, 'Children');

for i = 1:numel(h_patch)
    y = get(h_patch(i), 'YData');
    s = size(y);
    set(h_patch(i), 'ZData', zlevel*ones(s),'FaceColor',fcolor);
end

if (orgholdstate == 0), hold off; end;




function s = scaleme(x,minsize)
% scaleme returns a uniformly scaled array the same size as the input
% array where the maximum is 10 times the minimum 
maxsize = 10*minsize;
minx    = min(x);
maxx    = max(x);
slope   = (maxsize-minsize)/(maxx-minx);
b       = minsize - slope*minx;

s = x*slope + b;



function plot_3D(obsstruct, pstruct)

if (pstruct.clim(1) == pstruct.clim(2))
   % If all the observations have the same value, setting the
   % colorbar limits is a real pain. Fundamentally, I am 
   % forcing the plot symbols to be the lowest color of the
   % colormap and setting the colorbar to have some more
   % colors 'on top' - that are never used.
   cmap = colormap;
   h = plot3(obsstruct.lons, obsstruct.lats, obsstruct.z, 'bo');
   set(h,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
   set(gca,'Clim',[pstruct.clim(1) pstruct.clim(2)+1])
   set(gca,'XGrid','on','YGrid','on','ZGrid','on')

else
   scatter3(obsstruct.lons, obsstruct.lats, obsstruct.z, ...
         pstruct.scalearray, obsstruct.obs, 'd', 'filled');
end

clim = get(gca,'CLim');

axis(pstruct.axis)

title( {pstruct.str1, pstruct.str3, pstruct.str2}, 'Interpreter','none','FontSize',16);
xlabel('longitude')
ylabel('latitude')

if     (obsstruct.Ztyp(1) == -2) % VERTISUNDEF     = -2
   zlabel('unspecified')
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
set(gca,'CLim',clim)
h = colorbar;
set(get(h,'YLabel'),'String',pstruct.colorbarstring,'Interpreter','none')




function plot_2D(obsstruct, pstruct)

axis(pstruct.axis); hold on; worldmap('light');

if (pstruct.clim(1) == pstruct.clim(2))
   cmap = colormap;
   h = plot(obsstruct.lons, obsstruct.lats, 'bo');
   set(h,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
   set(gca,'Clim',[pstruct.clim(1) pstruct.clim(2)+1])
   set(gca,'XGrid','on','YGrid','on')

else

   scatter(obsstruct.lons, obsstruct.lats, ...
         pstruct.scalearray, obsstruct.obs, 'd', 'filled');
end

clim = get(gca,'CLim');

title( {pstruct.str1, pstruct.str3, pstruct.str2}, 'Interpreter','none','FontSize',16);
xlabel('longitude')
ylabel('latitude')

set(gca,'CLim',clim)
h = colorbar;
set(get(h,'YLabel'),'String',pstruct.colorbarstring,'Interpreter','none')
hold off
