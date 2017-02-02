%% rms_cross_time

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

prfname = 'Prior_Diag.nc';
pofname = 'Posterior_Diag.nc';
trfname = 'True_State.nc';

if (exist(trfname,'file') ~= 2)
   error('%s does not exist.',trfname)
end
if (exist(pofname,'file') ~= 2)
   error('%s does not exist.',pofname)
end
if (exist(prfname,'file') ~= 2)
   error('%s does not exist.',prfname)
end

if (~ nc_isvar(pofname,'state'))
    disp('This is an old-school routine ...')
    error('%s does not have a "state" variable',pofname)
end

tlon       = nc_varget(pofname, 'XLON_d01');
tlat       = nc_varget(pofname, 'XLAT_d01');
level      = nc_varget(pofname, 'level_d01');
copy       = nc_varget(pofname, 'copy');
true_times = nc_varget(pofname, 'time');

we             = size(      tlon, 2);
sn             = size(      tlat, 1);
bt             = size(     level, 1);
ens_size       = size(      copy, 1) - 2;
num_true_times = size(true_times, 1);

dx = 200.0;

mean_ind = get_copy_index(pofname,'ensemble mean');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

for field_num = 1:9

   start_var = 1;
   nx        = we + 1;
   ny        = sn;
   var_units = 'U (m/s)';
   maxlev    = bt;
   maxval    = 8.0;

   if field_num > 1
      start_var = start_var + bt*(we + 1)*sn ;
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
      maxval    = 700.0;
   end
   if field_num > 4
      start_var = start_var + (bt + 1)*we*sn;
      var_units = 'T (K)';
      maxlev    = bt;
      maxval    = 4.0;
   end
   if field_num > 5
      start_var = start_var + bt*we*sn;
      var_units = 'MU (Pa)';
      maxlev    = 1;
      maxval    = 700.0;
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

   x       = 1:nx;
   rmse    = zeros(nx,num_true_times);
   field1d = zeros(maxlev,1);
   f_size  = nx*ny*maxlev;
   end_var = start_var + f_size - 1;

   for itime = num_true_times-1:num_true_times

      % Extract field

      count = [    1  1 (end_var - start_var + 1)];
      start = [itime  1 start_var] -1;
      field_vec_truth     = nc_varget(trfname, 'state', start, count);

      start = [itime mean_ind start_var];
      field_vec_prior     = nc_varget(prfname, 'state', start, count);
      field_vec_posterior = nc_varget(pofname, 'state', start, count);

      field_vec = field_vec_posterior - field_vec_truth;
      field3d   = reshape(field_vec, [nx, ny, maxlev]);

      for ix = 1:nx
         field1d(:) = field3d(ix,sn/2,:);
         rmse(ix,itime) = sqrt((field1d'*field1d)/maxlev);
         x(ix) = (ix-1)*dx;
      end

   end

   subplot(3,3,field_num);
   plot(x,rmse)
   title(var_units)

end

xlabel('distance (km)')

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
