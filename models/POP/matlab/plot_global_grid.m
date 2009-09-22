function plot_global_grid(n)
% plot_global_grid 
% 
% plot_global_grid(1); % to plot regular grid output
% plot_global_grid(2); % to plot dipole x3 grid output
% plot_global_grid(3); % to plot regular grid with same grid as x3 in SH

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if n == 1
   str0 = 'regular grid';
   u_out = load('regular_grid_u_out');
   t_out = load('regular_grid_t_out');
   nx = 111;
   ny = 317;
elseif n == 2
   str0 = 'dipole grid';
   u_out = load('dipole_x3_u_out');
   t_out = load('dipole_x3_t_out');
   nx = 100;
   ny = 116;
elseif n == 3
   str0 = 'regular grid in NH, dipole in SH';
   u_out = load('regular_griddi_u_out');
   t_out = load('regular_griddi_t_out');
   nx = 100;
   ny = 97;
end

u = reshape(u_out(:,3),ny,nx);
t = reshape(t_out(:,3),ny,nx);

str1 = sprintf('min/max is %f %f',min(u(:)), max(u(:)));
str2 = sprintf('min/max is %f %f',min(t(:)), max(t(:)));

figure(1); clf
imagesc(u); colorbar;
title({'U Grid',str0,str1})

figure(2); clf
imagesc(t); colorbar;
title({'T Grid',str0,str2})

