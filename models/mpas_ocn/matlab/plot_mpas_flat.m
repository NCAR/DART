function h = plot_mpas_flat(fname,varname,level)
%% plot_mpas_flat
% 
% fname   = '../data/mpas_output.2010-10-23_03:00:00.nc';
% varname = 'theta';
% level   = 4;
% plot_mpas_flat(fname,varname,level)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%% Get the information about the polygons.

nEdgesOnCell   = nc_varget(fname,'nEdgesOnCell');
verticesOnCell = nc_varget(fname,'verticesOnCell');
lonCell        = nc_varget(fname,'lonCell') * 180.0/pi;
latCell        = nc_varget(fname,'latCell') * 180.0/pi;
lonVertex      = nc_varget(fname,'lonVertex') * 180.0/pi;
latVertex      = nc_varget(fname,'latVertex') * 180.0/pi;
nCells         = size(latCell,1);

[start, count] = find_hyperslab(fname,varname,level);
c = nc_varget(fname,varname,start,count);

minfld   = min(c(:));
maxfld   = max(c(:));
fldrange = maxfld - minfld;

%% Each Cell has some number of vertices

nrows = max(nEdgesOnCell);
xpoly = NaN(nrows,nCells);
ypoly = NaN(nrows,nCells);

for iCell=1:nCells

   n = nEdgesOnCell(iCell);

%  xpoly(1:n,iCell) = lonVertex(verticesOnCell(iCell,1:n));
%  ypoly(1:n,iCell) = latVertex(verticesOnCell(iCell,1:n));

%  inds = find(abs(xpoly(:,iCell) - xpoly(1,iCell)) > 180.0);

%  if (~isempty(inds))
%     for i = 2:n
%        if (xpoly(i,iCell) > xpoly(1,iCell)) 
%           xpoly(i,iCell) = xpoly(i,iCell) - 360.0;
%        else
%           xpoly(i,iCell) = xpoly(i,iCell) + 360.0;
%        end
%     end
%  end

   for i = 1:n
      xpoly(i,iCell) = lonVertex(verticesOnCell(iCell,i));
      ypoly(i,iCell) = latVertex(verticesOnCell(iCell,i));
      if (i > 1) 
         if (abs(xpoly(i,iCell) - xpoly(1,iCell)) > 180.0) 
            if (xpoly(i,iCell) > xpoly(1,iCell)) 
               xpoly(i,iCell) = xpoly(i,iCell) - 360.0;
            else
               xpoly(i,iCell) = xpoly(i,iCell) + 360.0;
            end
         end
      end
   end

   % patch the pentagons up to hexagons
   if (n < nrows) 
      xpoly(n+1:nrows,iCell) = xpoly(n,iCell);
      ypoly(n+1:nrows,iCell) = ypoly(n,iCell);
   end

end

figure(1)
h = patch(xpoly,ypoly,c);
set(h,'LineStyle','none');
colorbar
axis image

% figure(2)
% m_proj('set','Equidistant Cylindrical');
% h = m_patch(xpoly,ypoly,c);
% set(h,'LineStyle','none');
% colorbar


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
