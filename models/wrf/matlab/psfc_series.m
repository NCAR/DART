%% psfc_series

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

field_name = 'MU';

fname = 'psfc';

times    = nc_varget(fname,     'time'); Ntimes = size(times, 1);
xlon     = nc_varget(fname, 'XLON_d01');  we    = size( xlon, 2);
xlat     = nc_varget(fname, 'XLAT_d01');  sn    = size( xlat, 1);
level    = nc_varget(fname,'level_d01');  bt    = size(level, 1);
dt       = nc_varget(fname,       'DT');

f_size = we*sn;
ftime  = Ntimes - 1;
pane   = 1;
x      = 0:ftime-1;
rmse   = x;

for itime = 1:ftime
   if ~strcmp(times(itime,:),times(itime+1,:))

      %% Extract two adjacent timesteps of the field and manipulate

      start = [itime    1  1] -1;
      count = [    2   -1 -1];

      field12 = nc_varget(fname, field_name, start, count);
      field   = (field12(2,:,:) - field12(1,:,:))/dt;

      field_vec  = reshape(field,f_size,1);
      rmse(pane) = sqrt((field_vec'*field_vec)/(f_size));

      pane = pane + 1;

   end
end

pane = pane - 1;

x         = x*dt;
time_unit = 'seconds';

if (max(x) > 60)
   x = x/60;
   time_unit = 'minutes';
   if (max(x) > 60)
      x = x/60;
      time_unit = 'hours';
      if (max(x) > 24)
         x = x/24;
         time_unit = 'days';
      end
   end
end

%----------------------------------------------------------------------
figure(1); clf; orient landscape; wysiwyg;
%----------------------------------------------------------------------

axes('FontSize', 14)

plot(x(1:pane),rmse(1:pane),'k')
%plot(x,rmse,x,cv5_v1_lhalf,x,cv5_v1_l1,x,cv3)

%plot(x(1:pane),rmse_3dvar(1:pane),'b-','LineWidth',2)
%hold on
%plot(x(1:pane),rmse_enkf(1:pane),'r-','LineWidth',2)
%hold on
%plot(x(1:pane),rmse(1:pane),'k-','LineWidth',2)

%plot(x,rmse,'LineWidth',2)
xlabel(time_unit,'Fontsize',18)
ylabel('Pa/s','Fontsize',18)

%     legend('CV5 var-2.0 len-0.5','CV5 var-1.0 len-0.5','CV5 var-1.0 len-1.0','CV3')
%legend('3D-Var','EnKF','no assim')

mean(rmse(1:pane))
mean(rmse(1:72:pane))

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
