function plot_grid(fname)
%% plot_grid ... plots the ulat,ulon,ulev variables from a netcdf file.
% 
% fname = 'gcom_restart.nc';
% plot_grid(fname)

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

ulon3D = nc_varget(fname,'ulon');
ulat3D = nc_varget(fname,'ulat');
ulev3D = nc_varget(fname,'ulev');

deg2rad = pi/180.0;

[nlev, nlat, nlon] = size(ulon3D);

inds = find(ulon3D <= 0);
ulon3D(inds) = ulon3D(inds) + 360.0;

% Just pick off the 10th level.
Ulon = squeeze(ulon3D(10,:,:));
Ulat = squeeze(ulat3D(10,:,:));
Ulev = squeeze(ulev3D(10,:,:));

lon = Ulon(:);
lat = Ulat(:);
lev = Ulev(:);

figure(1); clf; orient landscape
   plot3(lon,lat,lev,'ko')
   title('U gridpoints for level 10')
   xlabel('longitude (degrees)')   % fixme ... read the attributes and plot
   ylabel('latitude (degrees)')    % fixme ... read the attributes and plot
   zlabel('depth (meters)')        % fixme ... read the attributes and plot

figure(2); clf; orient landscape

    iN = 64;
    jN = 16;
    kN = 10;

    fprintf('%20.15f %20.15f %20.15f\n', ulon3D(kN,jN,iN)        , ulat3D(kN,jN,iN)        , ulev3D(kN,jN,iN))
    fprintf('%20.15f %20.15f %20.15f\n', ulon3D(kN,jN,iN)*deg2rad, ulat3D(kN,jN,iN)*deg2rad, ulev3D(kN,jN,iN))

    my3Dplot(iN,jN,kN,ulon3D,ulat3D,ulev3D)

figure(3); clf; orient landscape

    my2Dplot(iN,jN,kN,ulon3D,ulat3D)



function my2Dplot(ilon,ilat,ilev,ulon,ulat)

   lon1 = ilon-3; lonN = ilon+3;
   lat1 = ilat-3; latN = ilat+3;

   lons = ulon(ilev, lat1:latN, lon1:lonN);
   lats = ulat(ilev, lat1:latN, lon1:lonN);

   h1 = plot(lons(:), lats(:), 'ko');
   set(h1,'MarkerSize',1)
   axis equal

   for i = lon1:lonN
   for j = lat1:latN
   
      str   = sprintf('''%d,%d''',i,j);
      mystr = sprintf('h1 = text(%f,%f,%s);',ulon(ilev,j,i),ulat(ilev,j,i),str);
      eval(mystr)
      set(h1,'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',10,'Color','k')
        
   end
   end

 
function my3Dplot(ilon,ilat,ilev,ulon,ulat,ulev)

   lon1 = ilon-3; lonN = ilon+3;
   lat1 = ilat-3; latN = ilat+3;
   lev1 = ilev-1; levN = ilev+1;

   lons = ulon(lev1:levN, lat1:latN, lon1:lonN);
   lats = ulat(lev1:levN, lat1:latN, lon1:lonN);
   levs = ulev(lev1:levN, lat1:latN, lon1:lonN);

   h1 = plot3(lons(:), lats(:), levs(:), 'ko');
   set(h1,'MarkerSize',10)

   hold on;

   h2 = plot3(ulon(ilev,ilat,ilon), ulat(ilev,ilat,ilon), ulev(ilev,ilat,ilon), 'ro');
   set(h2,'MarkerSize',10,'MarkerFaceColor','r')

   title({'GCOM U grid layout','64,16,10 in red'})
   xlabel('longitude (degrees)')
   ylabel('latitude (degrees)')
   zlabel('depth (m)')
   set(gca,'box','off')

   hold off;

% TJH FIXME ... it would be cool to have the same plot - in radians - with the 
% vert_normalization_height applied - and another with the anything within the 
% cutoff radius highlighted.

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

