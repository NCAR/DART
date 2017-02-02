function plot_global_grid(ugrid_file, tgrid_file)
%% plot_global_grid 
% 
% Example 1: plot the regular grid output
%
% ugrid_file = 'dipole_u_out_data';
% tgrid_file = 'dipole_t_out_data';
% plot_global_grid(ugrid_file, tgrid_file)
%
%
% plot_global_grid(1); % to plot regular grid output
% plot_global_grid(2); % to plot dipole x3 grid output
% plot_global_grid(3); % to plot regular grid with same grid as x3 in SH

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

u = read_file(ugrid_file);
t = read_file(tgrid_file);

function chunk = read_file(fname)

if (exist(fname,'file') ~=2 )
   error('%s does not exist',fname) 
end

chunk.fname  = fname;
mydata       = load(fname);
chunk.lons   = mydata(:,1);
chunk.lats   = mydata(:,2);
chunk.vals   = mydata(:,3);
chunk.nx     = sum(abs(diff(chunk.lats)) > 20) + 1;  % looking for big jumps
chunk.ny     = length(chunk.lons) / chunk.nx;
chunk.datmat = reshape(chunk.vals, chunk.ny, chunk.nx);
chunk.string = sprintf('min/max is %f %f',min(chunk.vals), max(chunk.vals));

figure;

% imagesc(chunk.lons, chunk.lats, chunk.datmat)  ! dipole grids have lon(1) ~= 0
imagesc(chunk.datmat) 
set(gca,'YDir','normal'); colorbar
h = title({chunk.fname, chunk.string}); set(h,'Interpreter','none')
xlabel(sprintf('%d longitudes',chunk.nx))
ylabel(sprintf('%d latitudes',chunk.ny))

figure;

h = plot(chunk.lons, chunk.lats,'.'); set(h,'MarkerSize',0.2)
h = title(chunk.fname); set(h,'Interpreter','none')
xlabel(sprintf('%d longitudes',chunk.nx))
ylabel(sprintf('%d latitudes',chunk.ny))

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
