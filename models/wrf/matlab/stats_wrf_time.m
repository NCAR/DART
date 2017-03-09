%% stats_wrf_time

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

[num_domains,~] = nc_dim_info(prfname,'domain');

if (num_domains > 1)

   disp(['Number of domains: ',int2str(num_domains)])
   id = input('Input domain id: ');

else

   id = 1;

end

[sls     , ~] = nc_dim_info(prfname,['soil_layers_stag_d0',int2str(id)]);
[we      , ~] = nc_dim_info(prfname,['west_east_d0',       int2str(id)]);
[sn      , ~] = nc_dim_info(prfname,['south_north_d0',     int2str(id)]);
[bt      , ~] = nc_dim_info(prfname,['bottom_top_d0',      int2str(id)]);
[ens_size, ~] = nc_dim_info(prfname,'member');

true_times     = nc_varget(trfname, 'time');
num_true_times = size(true_times, 1);
true_times     = true_times - true_times(1);
time_unit      = 'days';

if (true_times(num_true_times) < 1.0)
     true_times = true_times*24;
     time_unit  = 'hours';
end
if (true_times(num_true_times) < 1.0)
     true_times = true_times*60;
     time_unit  = 'minutes';
end
if (true_times(num_true_times) < 1.0)
     true_times = true_times*60;
     time_unit  = 'seconds';
end

true_ind = get_copy_index(trfname,'true state');
mean_ind = get_copy_index(prfname,'ensemble mean');
sprd_ind = get_copy_index(prfname,'ensemble spread');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.7*scrsz(3) 0.7*scrsz(4)])

axes('FontSize',12)

for field_num = 1:9

   nx        = we + 1;
   ny        = sn;
   var_units = 'U (m/s)';
   var_name  = ['U_d0',int2str(id)];
   maxlev    = bt;
   maxval    = 2.0;   % or 8.0

   if field_num > 1
      nx = we;
      ny        = sn;
      var_units = 'W (m/s)';
      var_name  = ['W_d0',int2str(id)];
      maxlev    = bt + 1;
      maxval    = 0.3; %maxval = 0.03;
   end
   if field_num > 2
      var_units = 'PH (m^2/s^2)';
      var_name  = ['PH_d0',int2str(id)];
      maxval    = 100.0; %maxval = 700.0;
   end
   if field_num > 3
      var_units = 'T (K)';
      var_name  = ['T_d0',int2str(id)];
      maxlev    = bt;
      maxval    = 1.0; %maxval = 4.0;
   end
   if field_num > 4
      var_units = 'MU (Pa)';
      var_name  = ['MU_d0',int2str(id)];
      maxlev    = 1;
      maxval    = 100.0; %maxval = 700.0;
   end
   if field_num > 5
      var_units = 'TSLB (K)';
      var_name  = ['TSLB_d0',int2str(id)];
      maxlev    = sls;
      maxval    = 0.1; %maxval = 4.0;
   end
   if field_num > 6
      var_units = 'QV (kg/kg)';
      var_name  = ['QVAPOR_d0',int2str(id)];
      maxlev    = bt;
      maxval    = 0.0005; %maxval = 0.001;
   end
   if field_num > 7
      var_units = 'QC (kg/kg)';
      var_name  = ['QCLOUD_d0',int2str(id)];
      maxval    = 0.00007;
   end
   if field_num > 8
      var_units = 'QR (kg/kg)';
      var_name  = ['QRAIN_d0',int2str(id)];
      maxval    = 0.0002; %maxval = 0.00007;
   end

   x      = 1:2*num_true_times;
   rmse   = x;
   spread = x;
   E2     = x-x;

   rms_mem = zeros(2*num_true_times,ens_size);

   f_size  = nx*ny*maxlev;

   for itime = 1:num_true_times

      if maxlev > 1
         start_t = [itime true_ind  1  1  1] - 1;
         start_m = [itime mean_ind  1  1  1] - 1;
         start_s = [itime sprd_ind  1  1  1] - 1;
         count   = [    1        1 -1 -1 -1];
      else
         start_t = [itime true_ind -1 -1] - 1;
         start_m = [itime mean_ind -1 -1] - 1;
         start_s = [itime sprd_ind -1 -1] - 1;
         count   = [    1        1 -1 -1];
   end

   %% Extract field

   field_vec_truth = reshape(nc_varget(trfname, var_name, start_t, count), f_size,1);
   field_vec_prior = reshape(nc_varget(prfname, var_name, start_m, count), f_size,1);
   field_vec       = reshape(nc_varget(prfname, var_name, start_s, count), f_size,1);

   spread(2*itime-1) = sqrt((field_vec'*field_vec)/(f_size));
   field_vec         = field_vec_prior - field_vec_truth;
   rmse(2*itime-1)   = sqrt((field_vec'*field_vec)/(f_size));

   x(2*itime-1)      = true_times(itime);

   field_vec_poste = reshape(nc_varget(pofname, var_name, start_m, count), f_size,1);
   field_vec       = reshape(nc_varget(pofname, var_name, start_s, count), f_size,1);

   spread(2*itime) = sqrt((field_vec'*field_vec)/(f_size));
   field_vec       = field_vec_poste - field_vec_truth;
   rmse(2*itime)   = sqrt((field_vec'*field_vec)/(f_size));

   x(2*itime) = true_times(itime);

   %% Loop over each member 

   for imem = 1:ens_size

      memstring = sprintf('ensemble member %d',imem);
      memindex  = get_copy_index(prfname,memstring);

      if maxlev > 1
         start = [itime memindex  1  1  1] - 1;
         count = [    1        1 -1 -1 -1];
      else
         start = [itime memindex  1  1] - 1;
         count = [    1        1 -1 -1];
      end

      field_vec_prior         = reshape(nc_varget(prfname,var_name,start,count),f_size,1);
      field_vec               = field_vec_prior - field_vec_truth;
      rms_mem(2*itime-1,imem) = sqrt((field_vec'*field_vec)/(f_size));

      field_vec_poste         = reshape(nc_varget(pofname,var_name,start,count),f_size,1);
      field_vec               = field_vec_poste - field_vec_truth;
      rms_mem(2*itime  ,imem) = sqrt((field_vec'*field_vec)/(f_size));

      E2(2*itime-1) = E2(2*itime-1) + rms_mem(2*itime-1,imem);
      E2(2*itime  ) = E2(2*itime  ) + rms_mem(2*itime  ,imem);
   
   end

end

   subplot(3,3,field_num);

   %if (ens_size > 0.0)
   %     E2 = E2/ens_size;
   %     plot(x,rmse,x,spread,x,E2,'--m')
   %else
        plot(x,rmse,x,spread,'LineWidth',2)
   %     plot(x,rmse,'--k',x,rmse_3dvar(:,field_num),':k',x,rmse_no_assim(:,field_num),'-.k','LineWidth',2)
   %     plot(x,rmse,x,rmse_OSSE31(:,field_num),x,rmse_OSSE30(:,field_num),x,rmse_OSSE26(:,field_num),x,rmse_OSSE32(:,field_num),x,rmse_OSSE33(:,field_num),x,rmse_OSSE34(:,field_num),'LineWidth',2)
   %end

   %rmse_OSSE34(:,field_num) = rmse;

        title(var_units,'Fontsize',12)

   %     xlabel(time_unit)

        room = (x(2*num_true_times)-x(1)+1)/10;

   %     axis ([(x(1)-room) (x(2*num_true_times)+room) min(min(rmse,spread)) max(max(rmse,spread))])
        axis ([(x(1)-room) (x(2*num_true_times)+room) 0.0 maxval])

   %     if (field_num == 6)
   %        legend('EnKF','3D-Var','no assim')
   %        legend('EnKF','3D-Var')
   %     end

end

%legend('Clim BC','20% 3D-Var BC','50% 3D-Var BC','100% 3D-Var BC','120% 3D-Var BC','150% 3D-Var BC','200% 3D-Var BC')

%if (ens_size > 0.0)
%     legend('RMS error','Spread','E2')
%else
     legend('RMS error','Spread')
%end
xlabel(time_unit,'Fontsize',12)

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
