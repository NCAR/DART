function plot_interp_diffs(n)
% plot_interp_diffs - script to examine the interpolation tests.
% 
% plot_interp_diffs(1) ;      % difference regular grids
% plot_interp_diffs(2) ;      % difference dipole x3 grids
% plot_interp_diffs(3) ;      % regular grids with same grid as x3 in SH

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

%n = input(['Input 1 to difference regular grids .\n' ...
%           'Input 2 to difference dipole x3 grids . \n' ...
%           'Input 3 to diffference regular grids with same grid as x3 in SH\n'])

if n == 1
   u_data = load('regular_grid_u_data');
   t_data = load('regular_grid_t_data');
   u_out  = load('regular_grid_u_out');
   t_out  = load('regular_grid_t_out');
   nx = 111;
   ny = 317;
elseif n == 2
   u_data = load('dipole_x3_u_data');
   t_data = load('dipole_x3_t_data');
   u_out  = load('dipole_x3_u_out');
   t_out  = load('dipole_x3_t_out');
   nx = 100;
   ny = 116;
elseif n == 3
   u_data = load('regular_griddi_u_data');
   t_data = load('regular_griddi_t_data');
   u_out  = load('regular_griddi_u_out');
   t_out  = load('regular_griddi_t_out');
   nx = 100;
   ny = 97;
else
   error('Sorry, choices are only 1, 2, or 3 ... you entered %d',n)
end

u  = reshape( u_out(:,3),ny,nx);
ub = reshape(u_data(:,3),ny,nx);
t  = reshape( t_out(:,3),ny,nx);
tb = reshape(t_data(:,3),ny,nx);

dif  = u - ub;
tdif = t - tb;

umin = min(dif(:));
umax = max(dif(:));
str1 = sprintf('min/max of difference is %f %f',umin,umax);

tmin = min(tdif(:));
tmax = max(tdif(:));
str2 = sprintf('min/max of difference is %f %f',tmin,tmax);

figure(1); clf;
imagesc(dif)
colorbar 
title({'U Difference',str1})

figure(2); clf;
imagesc(tdif)
colorbar
title({'T Difference',str2})
