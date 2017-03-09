%% correl

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear

fname = 'Prior_Diag.nc';

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

if (~ nc_isvar(fname,'state'))
    disp('This is an old-school routine ...')
    error('%s does not have a "state" variable',fname)
end

tlon  = ncread(fname,  'XLON_d01'); we = size( tlon, 2);
tlat  = ncread(fname,  'XLAT_d01'); sn = size( tlat, 1);
level = ncread(fname, 'level_d01'); bt = size(level, 1);

[ens_size,~] = nc_dim_info(fname,'member');

vartype = input('Input variable type for correlation, 1=U, 2=V, 3=W, 4=GZ, 5=T, 6=MU, 7=QV, 8=QC, 9=QR: ');
vartime = input('Input variable time (index): ');
vari    = input('Input variable x location : ');
varj    = input('Input variable y location : ');

if vartype == 6
   vark = 1;
else
   vark = input('Input variable z location : ');
end

start_var = 1;
nx        = we + 1;
ny        = sn;
if vartype > 1
   start_var = start_var + bt*(we + 1)*sn;
   nx        = we;
   ny        = sn + 1;
end
if vartype > 2
   start_var = start_var + bt*we*(sn + 1);
   nx        = we;
   ny        = sn;
end
if vartype > 3
   start_var = start_var + (bt + 1)*we*sn;
end
if vartype > 4
   start_var = start_var + (bt + 1)*we*sn;
end
if vartype > 5
   start_var = start_var + bt*we*sn;
end
if vartype > 6
   start_var = start_var + we*sn;
end
if vartype > 7
   start_var = start_var + bt*we*sn*(vartype-7);
end

start_var = start_var + vari + varj*nx + (vark-1)*nx*ny;
count     = [1 1 1];
ens_var  = zeros(1,ens_size);

for imem = 1:ens_size

   memstring = sprintf('ensemble member %d',imem);
   memindex  = get_member_index(fname,memstring);

   start = [vartime memindex start_var];
   ens_var(imem) = ncread(fname, 'state', start, count);
end

iso = 0.1:0.1:1;

%% Select field to plot (U, V, W, GZ, T, MU, QV, QC, QR)

field_num = input('Input field type, 1=U, 2=V, 3=W, 4=GZ, 5=T, 6=MU, 7=QV, 8=QC, 9=QR: ');

stime = input('Initial time : ');
ftime = input('End time : ');

%% Get level for free atmosphere fields
if field_num == 6
   field_level = 1;
else
   field_level = input('Input level: ');
end

start_var = 1;
nx        = we + 1;
ny        = sn;
var_units = 'U (m/s)';
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
end
if field_num > 3
   start_var = start_var + (bt + 1)*we*sn;
   var_units = 'GZ (m^2/s^2)';
end
if field_num > 4
   start_var = start_var + (bt + 1)*we*sn;
   var_units = 'T (K)';
end
if field_num > 5
   start_var = start_var + bt*we*sn;
   var_units = 'MU (Pa)';
end
if field_num > 6
   start_var = start_var + we*sn;
   var_units = 'QV (kg/kg)';
end
if field_num > 7
   start_var = start_var + bt*we*sn*(field_num-7);
   var_units = 'QC (kg/kg)';
end
if field_num > 8
   var_units = 'QR (kg/kg)';
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.9*scrsz(4) 0.9*scrsz(4)])

m         = ceil(sqrt(ftime-stime+1));
start_var = start_var + nx*ny*(field_level - 1);
end_var   = start_var + nx*ny - 1;
pane      = 1;

for itime = stime:ftime

   plot_title = [var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)];

   %% Extract field

   field_vec   = zeros(nx*ny,ens_size);
   correlation = zeros(nx*ny,1);
   count       = [1 1 1];

   for imem = 1:ens_size
      memstring         = sprintf('ensemble member %d',imem);
      memindex          = get_copy_index(fname,memstring);
      start             = [itime memindex start_var];
      state_vec_prior   = ncread(fname, 'state', start, count);
      field_vec(:,imem) = state_vec_prior;
   end

   iplot = 1;

   for i=1:nx*ny
      if ( (cov(ens_var) ~= 0.0) && (cov(field_vec(i,:)) ~= 0.0) )
         corrmat = corrcoef(ens_var,field_vec(i,:));
         correlation(i) = corrmat(1,2);
      else
         iplot=0;
      end
   end

   if iplot

      field = reshape(correlation, [nx, ny]);

      % Plot field

      subplot(m,m,pane);

      %nc=5

      %colormap = (prism(nc))

      if field_num > 2
         [C,h] = contour(tlon,tlat, field', iso );
         hold on
         [Cm,hm] = contour (tlon,tlat, field', -iso, ':');
         plot(tlon(vari+1),tlat(varj+1),'kx','LineWidth',2,'MarkerSize',15)
      else
         %[C, h] = contourf(field');
         [C,h] = contour (field', iso);
         hold on
         [Cm,hm] = contour (field', -iso, '--');
         plot(vari+1,varj+1,'kx','LineWidth',2,'MarkerSize',15)
      end
      title(plot_title)
      %colorbar('vert')
      clabel(C, h);
      clabel(Cm, hm);

   end

   pane = pane + 1;

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
