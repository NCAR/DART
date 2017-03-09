%% stats_wrf_time_vect_point

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

tlon  = nc_varget(prfname,  'XLON_d01'); we        = size(  tlon, 2);
tlat  = nc_varget(prfname,  'XLAT_d01'); sn        = size(  tlat, 1);
level = nc_varget(prfname, 'level_d01'); bt        = size( level, 1);
times = nc_varget(prfname,      'time'); num_times = size( times, 1);

[ens_size,~] = nc_dim_info(prfname,'member');

state_vec_truth = nc_varget(trfname, 'state');
state_vec_prior = nc_varget(prfname, 'state');
state_vec_poste = nc_varget(pofname, 'state');

dim = size(state_vec_truth);

len = dim(2);

state_vec_truth = reshape(state_vec_truth, [num_times, 1, len]);

true_ind = get_copy_index(prfname,'true state');
mean_ind = get_copy_index(prfname,'ensemble mean');
sprd_ind = get_copy_index(prfname,'ensemble spread');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

start_var = 1;
start_var = start_var + bt*(we + 1)*sn;
start_var = start_var + bt*we*(sn + 1);
%start_var = start_var + (bt + 1)*we*sn;
%start_var = start_var + (bt + 1)*we*sn;
%start_var = start_var + bt*we*sn;
%start_var = start_var + we*sn;
%start_var = start_var + bt*we*sn;
nx     = we;
ny     = sn;
maxlev = bt+1;
[maxqr,point] = max(state_vec_truth(num_times,1, start_var : start_var + nx*ny*maxlev - 1));

for field_num = 1:9

   start_var = 1;
   nx        = we + 1;
   ny        = sn;
   var_units = 'U (m/s)';
   maxlev    = bt;
   maxval    = 7.0;

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
      maxval    = 0.02;
   end
   if field_num > 3
      start_var = start_var + (bt + 1)*we*sn;
      var_units = 'GZ (m^2/s^2)';
      maxval    = 500.0;
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
      maxval    = 600.0;
   end
   if field_num > 6
      start_var = start_var + we*sn;
      var_units = 'QV (kg/kg)';
      maxlev    = bt;
      maxval    = 0.0007;
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

   x         = 1:2*num_times;
   rmse      = x;
   rmse_2    = x;
   spread    = x;
   spread_2  = x;
   start_var = start_var + point-1;
   rms_mem   = zeros(2*num_times,ens_size);

   for itime = 1:num_times

      %% Extract field

      field_vec = state_vec_prior(itime, mean_ind, start_var) - ...
                  state_vec_truth(itime, true_ind, start_var);

      rmse(  2*itime-1) = sqrt((field_vec*field_vec'));
      spread(2*itime-1) = state_vec_prior(itime, sprd_ind, start_var);
      x(     2*itime-1) = times(itime,1)*24; % time in hours

      field_vec = state_vec_poste(itime, mean_ind, start_var) - ...
                  state_vec_truth(itime, true_ind, start_var);

      rmse(  2*itime  ) = sqrt((field_vec*field_vec'));
      spread(2*itime  ) = state_vec_poste(itime, sprd_ind, start_var);
      x(     2*itime  ) = times(itime,1)*24; % time in hours

      mean_2   = state_vec_prior(itime, 1, start_var);
      mean_2   = mean_2 - mean_2;
      spread_f = mean_2;

      for imem = 1:ens_size

         memstring = sprintf('ensemble member %d',imem);
         memindex  = get_copy_index(prfname,memstring);

         mean_2 = mean_2 + state_vec_prior(itime, memindex, start_var);

         field_vec = state_vec_prior(itime, memindex, start_var) - ...
                     state_vec_truth(itime,        1, start_var);

         rms_mem(2*itime-1,imem) = sqrt((field_vec*field_vec'));
          
         field_vec = state_vec_poste(itime, memindex, start_var) - ...
                     state_vec_truth(    itime,        1, start_var);

         rms_mem(2*itime  ,imem) = sqrt((field_vec*field_vec'));

      end

      mean_2 = mean_2/ens_size;

      field_vec         = mean_2 - state_vec_truth(itime,1, start_var);
      rmse_2(2*itime-1) = sqrt((field_vec*field_vec'));
      rmse_2(2*itime  ) = sqrt((field_vec*field_vec'));

      for imem = 1:ens_size
         memstring = sprintf('ensemble member %d',imem);
         memindex  = get_copy_index(prfname,memstring);
         spread_f  = spread_f + (state_vec_prior(itime, memindex, start_var) - mean_2).^2;
      end

      spread_f            = sqrt(spread_f/(ens_size-1));
      spread_2(2*itime-1) = mean(spread_f);
      spread_2(2*itime  ) = mean(spread_f);

   end

   subplot(3,3,field_num);
   plot(x,rmse,x,spread,x,rms_mem,'-r',x,rmse_2,'--y',x,spread_2,'--m')

   title(var_units)

   room = (x(2*num_times)-x(1))/10;

   axis ([(x(1)-room) (x(2*num_times)+room) 0.0 maxval])

end

legend('RMS error','Spread')
xlabel('hours')

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
