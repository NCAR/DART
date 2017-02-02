function plot_mpas_marble(fname,varname,level)
%% plot_mpas_marble
% 
% fname   = '../data/mpas_output.2010-10-23_03:00:00.nc';
% varname = 'theta';
% level   = 4;
% plot_mpas_marble(fname,varname,level)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%% Get the information about the triangles.
x   = nc_varget(fname,'xCell');
y   = nc_varget(fname,'yCell');
z   = nc_varget(fname,'zCell');
tri = nc_varget(fname,'cellsOnVertex');

%% Get the information about the field to be plotted.

[start, count] = find_hyperslab(fname,varname,level);
c = nc_varget(fname,varname,start,count);
h = trisurf(tri,x,y,z,c);

set(h,'LineStyle','none')
axis equal
axis off
colorbar
% shading interp
h = title({sprintf('%s level %d',varname,level),fname});
set(h,'Interpreter','none');

function [start, count] = find_hyperslab(fname,varname,level);
% helper to retrieve just the hyperslab of interest
varinfo = nc_getvarinfo(fname,varname);
ndims   = length(varinfo.Dimension);
start   = zeros(1,ndims);
count   = zeros(1,ndims);

for i = 1:ndims
   diminfo = nc_getdiminfo(fname,varinfo.Dimension{i});
   dimname = diminfo.Name;
   switch lower(dimname)
      case {'nvertlevels','nvertlevelsp1','nsoillevels'}
         start(i) = level-1;
         count(i) = 1;
      otherwise
         count(i) = varinfo.Size(i);
   end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
