% Just making the coordinate variables in the shapes I believe we need them.

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%dimensions:
         nxp1 = 97 ;
         nyp1 = 33 ;
         nzp1 = 33 ;
%        time = UNLIMITED ; // (1 currently)
         nx = 96 ;
         ny = 32 ;
         nz = 32 ;
%variables:
%        double ulon(nz, ny, nxp1) ;
%                ulon:units = "kilometers_east" ;
%                ulon:long_name = "x-coordinate of the grid nodes" ;
%        double ulat(nz, ny, nxp1) ;
%                ulat:units = "kilometers_east" ;
%                ulat:long_name = "x-coordinate of the grid nodes" ;
%        double ulev(nz, ny, nxp1) ;
%                ulev:units = "kilometers_east" ;
%                ulev:long_name = "x-coordinate of the grid nodes" ;

ulon = zeros(nxp1,ny,nz);
ulat = zeros(nxp1,ny,nz);
ulev = zeros(nxp1,ny,nz);

 nccreate('grid.nc','ulon',...
          'Dimensions',{'nxp1' nxp1 'ny' ny 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','ulat',...
          'Dimensions',{'nxp1' nxp1 'ny' ny 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','ulev',...
          'Dimensions',{'nxp1' nxp1 'ny' ny 'nz' nz},...
          'Format','classic');

ncwrite('grid.nc','ulon', ulon);
ncwrite('grid.nc','ulat', ulat);
ncwrite('grid.nc','ulev', ulev);
ncwriteatt('grid.nc','ulon','long_name','longitude');
ncwriteatt('grid.nc','ulon','units','degrees_east');
ncwriteatt('grid.nc','ulon','axis','X');
ncwriteatt('grid.nc','ulat','long_name','latitude');
ncwriteatt('grid.nc','ulat','units','degrees_north');
ncwriteatt('grid.nc','ulat','axis','Y');
ncwriteatt('grid.nc','ulev','long_name','depth');
ncwriteatt('grid.nc','ulev','units','meters');
ncwriteatt('grid.nc','ulev','axis','Z');
ncwriteatt('grid.nc','ulev','positive','up');

ncwriteatt('grid.nc','ulon','comment','longitude of the u component of velocity');
ncwriteatt('grid.nc','ulat','comment','latitude of the u component of velocity');
ncwriteatt('grid.nc','ulev','comment','depth of the u component of velocity');

vlon = zeros(nx,nyp1,nz);
vlat = zeros(nx,nyp1,nz);
vlev = zeros(nx,nyp1,nz);

 nccreate('grid.nc','vlon',...
          'Dimensions',{'nx' nx 'nyp1' nyp1 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','vlat',...
          'Dimensions',{'nx' nx 'nyp1' nyp1 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','vlev',...
          'Dimensions',{'nx' nx 'nyp1' nyp1 'nz' nz},...
          'Format','classic');

ncwrite('grid.nc','vlon', vlon);
ncwrite('grid.nc','vlat', vlat);
ncwrite('grid.nc','vlev', vlev);
ncwriteatt('grid.nc','vlon','long_name','longitude');
ncwriteatt('grid.nc','vlon','units','degrees_east');
ncwriteatt('grid.nc','vlat','long_name','latitude');
ncwriteatt('grid.nc','vlat','units','degrees_north');
ncwriteatt('grid.nc','vlev','long_name','depth');
ncwriteatt('grid.nc','vlev','units','meters');
ncwriteatt('grid.nc','vlon','axis','X');
ncwriteatt('grid.nc','vlat','axis','Y');
ncwriteatt('grid.nc','vlev','axis','Z');
ncwriteatt('grid.nc','vlev','positive','up');

ncwriteatt('grid.nc','vlon','comment','longitude of the v component of velocity');
ncwriteatt('grid.nc','vlat','comment','latitude of the v component of velocity');
ncwriteatt('grid.nc','vlev','comment','depth of the v component of velocity');

wlon = zeros(nx,ny,nzp1);
wlat = zeros(nx,ny,nzp1);
wlev = zeros(nx,ny,nzp1);

 nccreate('grid.nc','wlon',...
          'Dimensions',{'nx' nx 'ny' ny 'nzp1' nzp1},...
          'Format','classic');

 nccreate('grid.nc','wlat',...
          'Dimensions',{'nx' nx 'ny' ny 'nzp1' nzp1},...
          'Format','classic');

 nccreate('grid.nc','wlev',...
          'Dimensions',{'nx' nx 'ny' ny 'nzp1' nzp1},...
          'Format','classic');

ncwrite('grid.nc','wlon', wlon);
ncwrite('grid.nc','wlat', wlat);
ncwrite('grid.nc','wlev', wlev);
ncwriteatt('grid.nc','wlon','long_name','longitude');
ncwriteatt('grid.nc','wlon','units','degrees_east');
ncwriteatt('grid.nc','wlat','long_name','latitude');
ncwriteatt('grid.nc','wlat','units','degrees_north');
ncwriteatt('grid.nc','wlev','long_name','depth');
ncwriteatt('grid.nc','wlev','units','meters');
ncwriteatt('grid.nc','wlon','axis','X');
ncwriteatt('grid.nc','wlat','axis','Y');
ncwriteatt('grid.nc','wlev','axis','Z');
ncwriteatt('grid.nc','wlev','positive','up');
ncwriteatt('grid.nc','wlon','comment','longitude of the w component of velocity');
ncwriteatt('grid.nc','wlat','comment','latitude of the w component of velocity');
ncwriteatt('grid.nc','wlev','comment','depth of the w component of velocity');

lon = zeros(nx,ny,nz);
lat = zeros(nx,ny,nz);
lev = zeros(nx,ny,nz);

 nccreate('grid.nc','lon',...
          'Dimensions',{'nx' nx 'ny' ny 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','lat',...
          'Dimensions',{'nx' nx 'ny' ny 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','lev',...
          'Dimensions',{'nx' nx 'ny' ny 'nz' nz},...
          'Format','classic');

 ncwrite('grid.nc','lon', lon);
 ncwrite('grid.nc','lat', lat);
 ncwrite('grid.nc','lev', lev);
ncwriteatt('grid.nc','lon','long_name','longitude');
ncwriteatt('grid.nc','lon','units','degrees_east');
ncwriteatt('grid.nc','lat','long_name','latitude');
ncwriteatt('grid.nc','lat','units','degrees_north');
ncwriteatt('grid.nc','lev','long_name','depth');
ncwriteatt('grid.nc','lev','units','meters');
ncwriteatt('grid.nc','lon','axis','X');
ncwriteatt('grid.nc','lat','axis','Y');
ncwriteatt('grid.nc','lev','axis','Z');
ncwriteatt('grid.nc','lev','positive','up');
ncwriteatt('grid.nc','lon','comment','longitude of the grid cell center');
ncwriteatt('grid.nc','lat','comment','latitude of the grid cell center');
ncwriteatt('grid.nc','lev','comment','depth of the grid cell center');

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

