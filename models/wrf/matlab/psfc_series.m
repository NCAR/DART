%% psfc_series

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

field_name = 'MU';

fname = 'psfc';

nc = netcdf( [fname,'.nc'] , 'nowrite' ) ;

we = size(nc('west_east'),1);
sn = size(nc('south_north'),1);
dt = nc.DT(:);

close(nc);

f_size = we*sn;

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);

ftime = num_true_times - 1;

pane = 1;

x = [0:ftime-1];
rmse = x;

for itime = 1:ftime

if ~strcmp(true_times(itime,:),true_times(itime+1,:))

% Extract field

   field1 = getnc(fname, field_name,[itime -1 -1],[itime -1 -1],[1 1 1]);
   field2 = getnc(fname, field_name,[itime+1 -1 -1],[itime+1 -1 -1],[1 1 1]);

   field = (field2 - field1)/dt;

   field_vec = reshape(field,f_size,1);

   rmse(pane) = sqrt((field_vec'*field_vec)/(f_size));

   pane = pane + 1;

end

end

pane = pane - 1;

x = x*dt;
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
