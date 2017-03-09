%% w2_wrf

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

tlon    = nc_varget(prfname,   'west_east'); we         = size(  tlon, 2);
tlat    = nc_varget(prfname, 'south_north'); sn         = size(  tlat, 1);
level   = nc_varget(prfname,  'bottom_top'); bt         = size( level, 1);
ttimes  = nc_varget(prfname,        'time'); num_ttimes = size(ttimes, 1);

[ens_size, ~] = nc_dim_info(prfname,'member');

true_ind = get_copy_index(prfname,'true state');
mean_ind = get_copy_index(prfname,'ensemble mean');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

field_num = 3;

start_var = 1;
start_var = start_var + bt*(we + 1)*sn ;
start_var = start_var + bt*we*(sn + 1);
nx        = we;
ny        = sn;
var_units = 'W (m/s)';
maxlev    = bt + 1;
maxval    = 0.03;

x       = 1:2*num_ttimes ;
estim   = x;
truth   = x;
rms_mem = zeros(2*num_ttimes ,ens_size);
f_size  = nx*ny*maxlev;
end_var = start_var + f_size - 1;

for itime = 1:num_ttimes 

   %% Extract fields

   start = [itime true_ind start_var] - 1;
   count = [itime        1 (end_var - start_var + 1)];

   field_vec_truth  = nc_varget(trfname, 'state', start, count);
   truth(2*itime-1) = sqrt((field_vec_truth'*field_vec_truth)/(f_size));
   truth(2*itime  ) = truth(2*itime-1);

   start = [itime mean_ind start_var] - 1;

   field_vec_prior  = nc_varget(prfname, 'state', start, count);
   estim(2*itime-1) = sqrt((field_vec_prior'*field_vec_prior)/(f_size));
   x(    2*itime-1) = ttimes (itime,1)*24; % time in hours

   field_vec_poste  = nc_varget(pofname, 'state', start, count);
   estim(2*itime)   = sqrt((field_vec_poste'*field_vec_poste)/(f_size));
   x(    2*itime)   = ttimes (itime,1)*24; % time in hours

   for imem = 1:ens_size

      start = [itime imem start_var] - 1;

      field_vec_prior         = nc_varget(prfname, 'state', start, count);
      field_vec               = field_vec_prior - field_vec_truth;
      rms_mem(2*itime-1,imem) = sqrt((field_vec'*field_vec)/(f_size));

      field_vec_poste       = nc_varget(pofname, 'state', start, count);
      field_vec             = field_vec_poste - field_vec_truth;
      rms_mem(2*itime,imem) = sqrt((field_vec'*field_vec)/(f_size));

      E2(2*itime-1) = E2(2*itime-1) + rms_mem(2*itime-1,imem);
      E2(2*itime  ) = E2(2*itime  ) + rms_mem(2*itime  ,imem);

   end

end


%subplot(3,3,field_num);

if (ens_size > 0.0)
     E2 = E2/ens_size;
     plot(x,truth,x,estim,x,E2,'--m')
else
     plot(x,truth,x,estim)
end

title(var_units)

room = (x(2*num_ttimes )-x(1))/10;

axis ([(x(1)-room) (x(2*num_ttimes )+room) 0.0 maxval])

legend('Truth','State Estimate')
xlabel('hours')

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
