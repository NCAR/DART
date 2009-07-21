clear; close all

n = input(['Input 1 to difference regular grids .\n' ...
          'Input 2 to difference dipole x3 grids . \n' ...
          'Input 3 to diffference regular grids with same grid as x3 in SH\n'])


if n == 1
   load regular_grid_u_data;
   u_data = regular_grid_u_data;
   load regular_grid_t_data;
   t_data = regular_grid_t_data;
   load regular_grid_u_out;
   u_out = regular_grid_u_out;
   load regular_grid_t_out;
   t_out = regular_grid_t_out;
   nx = 111;
   ny = 317;
elseif n == 2
   load dipole_x3_u_data;
   u_data = dipole_x3_u_data;
   load dipole_x3_t_data;
   t_data = dipole_x3_t_data;
   load dipole_x3_u_out;
   u_out = dipole_x3_u_out;
   load dipole_x3_t_out;
   t_out = dipole_x3_t_out;
   nx = 100;
   ny = 116;
elseif n == 3
   load regular_griddi_u_data;
   u_data = regular_griddi_u_data;
   load regular_griddi_t_data;
   t_data = regular_griddi_t_data;
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
      ub(i, j) = u_data(index, 3);
      dif(i, j) = u(i, j) - ub(i, j);
      t(i, j) = t_out(index, 3);
      tb(i, j) = t_data(index, 3);
      tdif(i, j) = t(i, j) - tb(i, j);
   end
end 

figure(1)
contour(dif', 20)
colorbar 
title 'U Difference'

figure(2)
contour(tdif', 20)
colorbar
title 'T Difference'

 


