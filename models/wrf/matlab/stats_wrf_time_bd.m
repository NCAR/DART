%% stats_wrf_time_bd

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

trfname = 'True_State.nc';
prfname = 'Prior_Diag.nc';
pofname = 'Posterior_Diag.nc';

if (exist(trfname,'file') ~= 2)
   error('%s does not exist.',trfname)
end
if (exist(prfname,'file') ~= 2)
   error('%s does not exist.',prfname)
end
if (exist(pofname,'file') ~= 2)
   error('%s does not exist.',pofname)
end

if (~ nc_isvar(prfname,'state'))
    disp('This is an old-school routine ...')
    error('%s does not have a "state" variable',prfname)
end

tlon  = nc_varget(prfname,  'XLON_d01'); we        = size( tlon, 2);
tlat  = nc_varget(prfname,  'XLAT_d01'); sn        = size( tlat, 1);
level = nc_varget(prfname, 'level_d01'); bt        = size(level, 1);
times = nc_varget(prfname,      'time'); num_times = size(times, 1);

[ens_size,~] = nc_dim_info(prfname,'member');
true_ind = get_copy_index(trfname,'true state');
mean_ind = get_copy_index(prfname,'ensemble mean');
sprd_ind = get_copy_index(prfname,'ensemble spread');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

for field_num = 1:9

   start_var = 1;
   nx        = we + 1;
   ny        = sn;
   var_units = 'U (m/s)';
   maxlev    = bt;
   maxval    = 6.0;

   if field_num > 1
      start_var = start_var + bt*(we + 1)*sn;
      nx        = we;
      ny        = sn + 1;
      var_units = 'V (m/s)';
   end
   if field_num > 2
      start_var = start_var + bt*we*(sn + 1);
      nx        = we;
      ny        = sn;
      var_units = 'W (m/s)';
      maxlev    = bt + 1;
      maxval    = 0.03;
   end
   if field_num > 3
      start_var = start_var + (bt + 1)*we*sn;
      var_units = 'GZ (m^2/s^2)';
      maxval    = 400.0;
   end
   if field_num > 4
      start_var = start_var + (bt + 1)*we*sn;
      var_units = 'T (K)';
      maxlev    = bt;
      maxval    = 3.0;
   end
   if field_num > 5
      start_var = start_var + bt*we*sn;
      var_units = 'MU (Pa)';
      maxlev    = 1;
      maxval    = 500.0;
   end
   if field_num > 6
      start_var = start_var + we*sn;
      var_units = 'QV (kg/kg)';
      maxlev    = bt;
      maxval    = 0.001;
   end
   if field_num > 7
      start_var = start_var + bt*we*sn*(field_num-7);
      var_units = 'QC (kg/kg)';
      maxval    = 0.00003;
   end
   if field_num > 8
      var_units = 'QR (kg/kg)';
      maxval    = 0.00007;
   end

   x      = 1:2*num_times;
   rmse   = x;
   spread = x;
   E2     = x-x;

   rms_mem = zeros(2*num_times,ens_size);
   f_size  = nx*ny*maxlev;
   end_var = start_var + f_size - 1;

   for itime = 1:num_times

   %% Extract field

   start = [itime true_ind start_var] - 1;
   count = [    1        1   (end_var - start_var + 1)];
   field_vec_truth = reshape(nc_varget(trfname,'state',start, count), [nx, ny, maxlev]);

   start = [itime mean_ind start_var] - 1;
   field_vec_prior = reshape(nc_varget(prfname,'state',start, count), [nx, ny, maxlev]);

   start = [itime sprd_ind start_var] - 1;
   field_vec       = reshape(nc_varget(prfname,'state',start, count), [nx, ny, maxlev]);

   spread(2*itime-1) = 0.0;

   for iz = 1:maxlev
      for ix = 1:nx
      spread(2*itime-1) = spread(2*itime-1) + field_vec(ix, 1,iz)^2.0;
      spread(2*itime-1) = spread(2*itime-1) + field_vec(ix,ny,iz)^2.0;
      end

      for iy = 2:ny-1
      spread(2*itime-1) = spread(2*itime-1) + field_vec( 1,iy,iz)^2.0;
      spread(2*itime-1) = spread(2*itime-1) + field_vec(nx,iy,iz)^2.0;
      end
   end

   spread(2*itime-1) = sqrt(spread(2*itime-1)/(maxlev*(2*nx+2*(ny-2))));

   field_vec = field_vec_prior - field_vec_truth;
   rmse(2*itime-1) = 0.0;

   for iz = 1:maxlev
      for ix = 1:nx
      rmse(2*itime-1) = rmse(2*itime-1) + field_vec(ix,1,iz)^2.0;
      rmse(2*itime-1) = rmse(2*itime-1) + field_vec(ix,ny,iz)^2.0;
      end

      for iy = 2:ny-1
      rmse(2*itime-1) = rmse(2*itime-1) + field_vec(1,iy,iz)^2.0;
      rmse(2*itime-1) = rmse(2*itime-1) + field_vec(nx,iy,iz)^2.0;
      end
   end

   rmse(2*itime-1) = sqrt(rmse(2*itime-1)/(maxlev*(2*nx+2*(ny-2))));

   x(2*itime-1) = times(itime,1)*24; % time in hours

   %% now for the posterior
   start = [itime mean_ind start_var] - 1;
   field_vec_poste = reshape(nc_varget(pofname, 'state', start, count), [nx, ny, maxlev]);

   start = [itime sprd_ind start_var] - 1;
   field_vec       = reshape(nc_varget(pofname, 'state', start, count), [nx, ny, maxlev]);

   spread(2*itime) = 0.0;

   for iz = 1:maxlev
      for ix = 1:nx
      spread(2*itime) = spread(2*itime) + field_vec(ix, 1,iz)^2.0;
      spread(2*itime) = spread(2*itime) + field_vec(ix,ny,iz)^2.0;
      end

      for iy = 2:ny-1
      spread(2*itime) = spread(2*itime) + field_vec( 1,iy,iz)^2.0;
      spread(2*itime) = spread(2*itime) + field_vec(nx,iy,iz)^2.0;
      end
   end

   spread(2*itime) = sqrt(spread(2*itime)/(maxlev*(2*nx+2*(ny-2))));

   field_vec = field_vec_poste - field_vec_truth;
   rmse(2*itime) = 0.0;

   for iz = 1:maxlev
      for ix = 1:nx
      rmse(2*itime) = rmse(2*itime) + field_vec(ix, 1,iz)^2.0;
      rmse(2*itime) = rmse(2*itime) + field_vec(ix,ny,iz)^2.0;
      end

      for iy = 2:ny-1
      rmse(2*itime) = rmse(2*itime) + field_vec( 1,iy,iz)^2.0;
      rmse(2*itime) = rmse(2*itime) + field_vec(nx,iy,iz)^2.0;
      end
   end

   rmse(2*itime) = sqrt(rmse(2*itime)/(maxlev*(2*nx+2*(ny-2))));

   x(2*itime) = times(itime,1)*24; % time in hours

   for imem = 1:ens_size

      memstring = sprintf('ensemble member %d',imem);
      memindex  = get_copy_index(fname, memstring);

      start = [itime memindex start_var] - 1;
      count = [itime memindex (end_var - start_var + 1)];

      field_vec_prior = nc_varget(prfname, 'state', start, count);
      field_vec       = field_vec_prior - field_vec_truth;
      rms_mem(2*itime-1,imem) = sqrt((field_vec'*field_vec)/(f_size));

      field_vec_poste = nc_varget(pofname, 'state', start, count);
      field_vec       = field_vec_poste - field_vec_truth;
      rms_mem(2*itime,imem) = sqrt((field_vec'*field_vec)/(f_size));

      E2(2*itime-1) = E2(2*itime-1) + rms_mem(2*itime-1,imem);
      E2(2*itime  ) = E2(2*itime  ) + rms_mem(2*itime  ,imem);

   end

end


subplot(3,3,field_num);

if (ens_size > 0.0)
     E2 = E2/ens_size;
     plot(x,rmse,x,spread,x,E2,'--m')
else
     plot(x,rmse,x,spread)
end

     title(var_units)

%     xlabel('hours')

     room = (x(2*num_times)-x(1))/10;

     axis ([(x(1)-room) (x(2*num_times)+room) 0.0 maxval])

end

legend('BD RMS error','BD Spread')
xlabel('hours')

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
