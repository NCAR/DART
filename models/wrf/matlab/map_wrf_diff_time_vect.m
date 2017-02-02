%% map_wrf_diff_time_vect

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

% Select field to plot (U, V, W, GZ, T, MU, QV, QC, QR)

field_num = input('Input field type, 1=U, 2=V, 3=W, 4=GZ, 5=T, 6=MU, 7=QV, 8=QC, 9=QR: ');

member = input('Input ensemble member index: ');

tlon  = nc_varget(prfname,  'XLON_d01'); we = size( tlon, 2);
tlat  = nc_varget(prfname,  'XLAT_d01'); sn = size( tlat, 1);
level = nc_varget(prfname, 'level_d01'); bt = size(level, 1);

stime = input('Initial time (index): ');
ftime = input('End time (index): ');

% Get level for free atmosphere fields
if field_num == 6
   field_level = 1;
else
   field_level = input('Input level: ');
end

start_var = 1;
nx        = we + 1;
ny        = sn;
var_units = 'U (m/s)';
iso       = 0.5:1:5;

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
   iso       = 0.01:0.01:0.1;
end
if field_num > 3
   start_var = start_var + (bt + 1)*we*sn;
   var_units = 'GZ (m^2/s^2)';
   iso       = 50:50:300;
end
if field_num > 4
   start_var = start_var + (bt + 1)*we*sn;
   var_units = 'T (K)';
   iso       = 0.5:0.5:5;
end
if field_num > 5
   start_var = start_var + bt*we*sn;
   var_units = 'MU (Pa)';
   iso       = 100:100:600;
end
if field_num > 6
   start_var = start_var + we*sn;
   var_units = 'QV (kg/kg)';
   iso       = 0.0001:0.0001:0.001;
end
if field_num > 7
   start_var = start_var + bt*we*sn*(field_num-7);
   var_units = 'QC (kg/kg)';
   iso       = 0.00001:0.00001:0.0001;
end
if field_num > 8
   var_units = 'QR (kg/kg)';
   iso       = 0.00001:0.00001:0.0001;
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.9*scrsz(4) 0.9*scrsz(4)])

m = ceil(sqrt(ftime-stime+1));

start_var = start_var + nx*ny*(field_level - 1);
end_var   = start_var + nx*ny - 1;

pane = 1;

for itime = stime:ftime

   plot_title = [var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)];

   %% Extract field
   count = [    1 1 (end_var-start_var+1)];
   tstart = [itime true_ind start_var] - 1;
   pstart = [itime   member start_var] - 1;

   state_vec_truth = nc_varget(trfname, 'state', tstart, count);
   state_vec_prior = nc_varget(prfname, 'state', pstart, count);
   state_vec_poste = nc_varget(pofname, 'state', pstart, count);

   %field_vec = state_vec_prior - state_vec_truth;
   field_vec = state_vec_poste - state_vec_prior;
   %field_vec = state_vec_poste;

   field = reshape(field_vec, [nx, ny]);

   % Plot field

   subplot(m,m,pane);

   %nc=5

   %colormap = (prism(nc))
   if field_num > 2
     [C,h] = contour(tlon,tlat, field', iso );
     hold on
     [Cm,hm] = contour (tlon,tlat, field', -iso, ':');
   else
     %[C, h] = contourf(field');
     [C,h] = contour (field', iso);
     hold on
     [Cm,hm] = contour (field', -iso, '--');
   end
   title(plot_title)
   %colorbar('vert')
   clabel(C, h);
   clabel(Cm, hm);

   pane = pane + 1;

end

% Loop for another try
%map_wrf;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
