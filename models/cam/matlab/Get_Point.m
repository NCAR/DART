function bob = Get_Point(filename,basevar,timeind,lat,lon,level,enssize)
%
% filename = '/project/dart/raeder/J/T85_3/01_04/Prior_Diag.nc';
% basevar = 'T';
% timeind = 2;
% level   = 7;
% lat = 40;      % index number, not direct lat value
% lon = 140;     % index number, not direct lon value
% enssize = 20;
% x = Get_point(filename,basevar,timeind,lat,lon,level,enssize);

% use some of the infinite numbers of options on getnc() to get only a
% hyperslab of the data in the first place rather than squeeze afterwards.

% assumes incoming data dims are: [timestep, ens_num, lat, lon, lev];

% first 2 copies are data mean/var; last 2 are inflation mean/var
bl_corner = [timeind, 3,         lat, lon, level];
ur_corner = [timeind, enssize+2, lat, lon, level];
squeeze_it = 1;

bob = getnc(filename, basevar, bl_corner, ur_corner, -1,-1,-1,-1, squeeze_it);

% doc for getnc()
%   function values = getnc(file, varid, corner, end_point, stride, order, ...
%                           change_miss, new_miss, squeeze_it, rescale_opts)
 

