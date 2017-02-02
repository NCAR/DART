%% stats_wrf_prof_mem

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

state_vec_truth = nc_varget(trfname, 'state');
state_vec_prior = nc_varget(prfname, 'state');
state_vec_poste = nc_varget(pofname, 'state');

dim = size(state_vec_truth);

len = dim(2);

state_vec_truth = reshape(state_vec_truth, [num_times, 1, len]);

member = input('Input ensemble member: ');
itime  = input('Time (index): ');

poste_tru = state_vec_poste(itime,member,:) - state_vec_truth(itime,1,:);
prior_tru = state_vec_prior(itime,member,:) - state_vec_truth(itime,1,:);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

for field_num = [1:5 7:9]

   start_var = 1;
   nx        = we + 1;
   ny        = sn;
   var_units = 'U (m/s)';
   maxlev    = bt;

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
   end
   if field_num > 3
      start_var = start_var + (bt + 1)*we*sn;
      var_units = 'GZ (m^2/s^2)';
   end
   if field_num > 4
      start_var = start_var + (bt + 1)*we*sn;
      var_units = 'T (K)';
      maxlev    = bt;
   end
   if field_num > 5
      start_var = start_var + bt*we*sn;
      var_units = 'MU (Pa)';
      maxlev    = 1;
   end
   if field_num > 6
      start_var = start_var + we*sn;
      var_units = 'QV (kg/kg)';
      maxlev    = bt;
   end
   if field_num > 7
      start_var = start_var + bt*we*sn*(field_num-7);
      var_units = 'QC (kg/kg)';
   end
   if field_num > 8
      var_units = 'QR (kg/kg)';
   end

   y          = 1:maxlev;
   rmse_prior = y;
   rmse_poste = y;

   for field_level = 1:maxlev

      % Extract field

      field_vec = prior_tru(start_var : start_var + nx*ny - 1);
      rmse_prior(field_level) = sqrt((field_vec*field_vec')/(nx*ny));

      field_vec = poste_tru(start_var : start_var + nx*ny - 1);
      rmse_poste(field_level) = sqrt((field_vec*field_vec')/(nx*ny));

      start_var = start_var + nx*ny;

   end

   subplot(3,3,field_num);
   plot(rmse_prior,y,rmse_poste,y)
   title(var_units)

end

legend('Prior','Posterior');

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
