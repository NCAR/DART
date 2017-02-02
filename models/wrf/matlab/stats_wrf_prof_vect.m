%% stats_wrf_prof_vect

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

tlon  = nc_varget(prfname,  'XLON_d01'); we = size( tlon, 2);
tlat  = nc_varget(prfname,  'XLAT_d01'); sn = size( tlat, 1);
level = nc_varget(prfname, 'level_d01'); bt = size(level, 1);

true_ind = get_copy_index(trfname,'true state');
mean_ind = get_copy_index(prfname,'ensemble mean');
sprd_ind = get_copy_index(prfname,'ensemble spread');

itime    = input('Time (index): ');

figure_title = ['Time: ',num2str(itime)];

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 0.9*scrsz(4)],'Name',figure_title)

for field_num = [1:5 7:9]

   start_var = 1;
   nx        = we + 1;
   ny        = sn;
   var_units = 'U (m/s)';
   maxlev    = bt;
   maxval    = 8.0;

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
      maxval    = 0.00004;
   end
   if field_num > 8
      var_units = 'QR (kg/kg)';
      maxval    = 0.00007;
   end

   f_size     = nx*ny;
   y          = 1:maxlev;
   rmse_prior = y;
   rmse_poste = y;
   sprd_prior = y;
   sprd_poste = y;

   for field_level = 1:maxlev

      %% Extract field

      start = [itime true_ind start_var] - 1;
      count = [    1        1     nx*ny];

      state_vec_truth = nc_varget(trfname, 'state', start, count);

      start = [itime mean_ind start_var] - 1;

      field_vec_prior = nc_varget(prfname, 'state', start, count);
      field_vec_poste = nc_varget(pofname, 'state', start, count);

      start = [itime sprd_ind start_var] - 1;

      field_vec               = nc_varget(prfname, 'state',start, count);
      sprd_prior(field_level) = sqrt((field_vec'*field_vec)/(f_size));
      field_vec               = nc_varget(pofname, 'state', start, count);
      sprd_poste(field_level) = sqrt((field_vec'*field_vec)/(f_size));

      field_vec               = field_vec_prior - state_vec_truth;
      rmse_prior(field_level) = sqrt((field_vec'*field_vec)/(f_size));
      field_vec               = field_vec_poste - state_vec_truth;
      rmse_poste(field_level) = sqrt((field_vec'*field_vec)/(f_size));

      start_var = start_var + nx*ny;

   end

   subplot(3,3,field_num);
   plot(rmse_prior,y,'b',rmse_poste,y,'g',sprd_prior,y,'--b',sprd_poste,y,'--g');
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
