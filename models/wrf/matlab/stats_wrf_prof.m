%% stats_wrf_prof

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

tlon  = nc_varget(prfname,  'XLON_d01'); we = size( tlon, 2);
tlat  = nc_varget(prfname,  'XLAT_d01'); sn = size( tlat, 1);
level = nc_varget(prfname, 'level_d01'); bt = size(level, 1);

true_ind = get_copy_index(trfname,'true state');
mean_ind = get_copy_index(prfname,'ensemble mean');
sprd_ind = get_copy_index(prfname,'ensemble spread');

itime = input('Time: ');

figure_title = ['Time: ',num2str(itime)];

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 0.9*scrsz(4)],'Name',figure_title)

for field_num = [1:5 7:9]

   nx        = we + 1;
   ny        = sn;
   var_units = 'U (m/s)';
   var_name  = 'U';
   maxlev    = bt;
   maxval    = 8.0;

   if field_num > 1
      nx = we;
      ny = sn + 1;
      var_units = 'V (m/s)';
      var_name = 'V';
   end
   if field_num > 2
      nx        = we;
      ny        = sn;
      var_units = 'W (m/s)';
      var_name  = 'W';
      maxlev    = bt + 1;
      maxval    = 0.03;
   end
   if field_num > 3
      var_units = 'GZ (m^2/s^2)';
      var_name  = 'PH';
      maxval    = 700.0;
   end
   if field_num > 4
      var_units = 'T (K)';
      var_name  = 'T';
      maxlev    = bt;
      maxval    = 4.0;
   end
   if field_num > 5
      var_units = 'MU (Pa)';
      var_name  = 'MU';
      maxlev    = 1;
      maxval    = 700.0;
   end
   if field_num > 6
      var_units = 'QV (kg/kg)';
      var_name  = 'QVAPOR';
      maxlev    = bt;
      maxval    = 0.001;
   end
   if field_num > 7
      var_units = 'QC (kg/kg)';
      var_name  = 'QCLOUD';
      maxval    = 0.00004;
   end
   if field_num > 8
      var_units = 'QR (kg/kg)';
      var_name  = 'QRAIN';
      maxval    = 0.00007;
   end

   f_size       = nx*ny;
   y            = 1:maxlev;
   rmse_prior   = y;
   rmse_poste   = y;
   spread_prior = y;
   spread_poste = y;

   for field_level = 1:maxlev

      count  = [1 1 1 -1 -1]; % time, copy, level, lat, lon

      %% Extract field

      start = [itime true_ind field_level  1  1] - 1;
      state_vec_truth = reshape(nc_varget(trfname, var_name, start, count),f_size,1);
   
      start = [itime mean_ind field_level  1  1] - 1;
      field_vec_prior = reshape(nc_varget(prfname, var_name, start, count),f_size,1);
      field_vec_poste = reshape(nc_varget(pofname, var_name, start, count),f_size,1);

      start = [itime sprd_ind field_level  1  1] - 1;
      bob = reshape(nc_varget(prfname, var_name, start, count),f_size,1);
      spread_prior(field_level) = sqrt((field_vec'*field_vec)/(f_size));

      bob = reshape(nc_varget(pofname, var_name, start, corner),f_size,1);
      spread_poste(field_level) = sqrt((field_vec'*field_vec)/(f_size));

      field_vec                 = field_vec_prior - state_vec_truth;
      rmse_prior(field_level)   = sqrt((field_vec'*field_vec)/(f_size));
      field_vec                 = field_vec_poste - state_vec_truth;
      rmse_poste(field_level)   = sqrt((field_vec'*field_vec)/(f_size));

   end

   subplot(3,3,field_num);
   plot(rmse_prior,y,'b',rmse_poste,y,'g',spread_prior,y,'--b',spread_poste,y,'--g');
   title(var_units);
   axis ([0.0 maxval 0.0 30.0])

end

legend('Prior','Posterior');

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
