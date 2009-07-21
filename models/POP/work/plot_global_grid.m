clear; close all

n = input(['Input 1 to plot regular grid output.\n' ...
          'Input 2 to plot dipole x3 grid output. \n' ...
          'Input 3 to plot regular grid with same grid as x3 in SH\n'])


if n == 1
   load regular_grid_u_out;
   u_out = regular_grid_u_out;
   load regular_grid_t_out;
   t_out = regular_grid_t_out;
   nx = 111;
   ny = 317;
elseif n == 2
   load dipole_x3_u_out;
   u_out = dipole_x3_u_out;
   load dipole_x3_t_out;
   t_out = dipole_x3_t_out;
   nx = 100;
   ny = 116;
elseif n == 3
   load regular_griddi_u_out;
   u_out = regular_griddi_u_out;
   load regular_griddi_t_out;
   t_out = regular_griddi_t_out;
   nx = 100;
   ny = 97;
end

index = 0;
for i = 1:nx
   for j = 1:ny
      index = index + 1;
      u(i, j) = u_out(index, 3);
      t(i, j) = t_out(index, 3);
   end
end 

figure(1)
contour(u')
title 'U Grid'

figure(2)
contour(t')
title 'T Grid'
      


