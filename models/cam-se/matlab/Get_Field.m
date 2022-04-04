function bob = Get_Field(filename,basevar,timeind,level,enssize)
%% x = Get_Field(filename,basevar,timeind,level);
% filename = '/project/dart/raeder/J/T85_3/01_04/Prior_Diag.nc';
% basevar = 'T';
% timeind = 2;
% level   = 7;

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% use some of the infinite numbers of options on getnc() to get only a
% hyperslab of the data in the first place rather than squeeze afterwards.

% assumes incoming data dims are: [timestep, ens_num, lat, lon, lev];

% first 2 copies are mean/var; last 2 are inflation mean/var
bl_corner = [timeind, 3,         -1, -1, level];
ur_corner = [timeind, enssize+2, -1, -1, level];
squeeze_it = 1;

bob = getnc(filename, basevar, bl_corner, ur_corner, -1,-1,-1,-1, squeeze_it);

% doc for getnc()
%   function values = getnc(file, varid, corner, end_point, stride, order, ...
%                           change_miss, new_miss, squeeze_it, rescale_opts)
 
% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
