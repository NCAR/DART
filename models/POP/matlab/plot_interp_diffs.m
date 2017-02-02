function plot_interp_diffs(ugrid_file, tgrid_file)
%% plot_interp_diffs - script to examine the interpolation tests.
% 
% Example 1: plot the regular grid output
%--------------------------------------------------
%    ugrid_file = 'regular_grid_u_data';
%    tgrid_file = 'regular_grid_t_data';
%    plot_interp_diffs(ugrid_file, tgrid_file)
%
% Example 2: plot the dipole grid output
%--------------------------------------------------
%    ugrid_file = 'dipole_grid_u_data';
%    tgrid_file = 'dipole_grid_t_data';
%    plot_interp_diffs(ugrid_file, tgrid_file)
%
% Example 3: plot the dipole grid output
%--------------------------------------------------
%    ugrid_file = 'regular_griddi_u_data';
%    tgrid_file = 'regular_griddi_t_data';
%    plot_interp_diffs(ugrid_file, tgrid_file)
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

u_org = read_file(ugrid_file);
t_org = read_file(tgrid_file);
u_new = read_file(sprintf('%s.out',ugrid_file));
t_new = read_file(sprintf('%s.out',tgrid_file));

udif = u_org.datmat - u_new.datmat;
tdif = t_org.datmat - t_new.datmat;

umin = min(udif(:));
umax = max(udif(:));
str1 = sprintf('min/max of difference is %f %f',umin,umax);

tmin = min(tdif(:));
tmax = max(tdif(:));
str2 = sprintf('min/max of difference is %f %f',tmin,tmax);

figure(1); clf;
imagesc(udif); set(gca,'YDir','normal')
colorbar 
title({'U Difference',str1})

figure(2); clf;
imagesc(tdif); set(gca,'YDir','normal')
colorbar
title({'T Difference',str2})


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

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
