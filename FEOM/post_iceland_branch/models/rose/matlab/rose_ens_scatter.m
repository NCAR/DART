% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
fname = input('Input file name for truth');
%fname = '../work/True_State.nc';
lon = getnc(fname, 'lon');
num_lon = size(lon, 1);
lat = getnc(fname, 'lat');
num_lat = size(lat, 1);
level = getnc(fname, 'lev');
num_level = size(level, 1);
times = getnc(fname, 'time');
num_times = size(times, 1);


%Load the ensemble file%
ens_fname = input('Input file name for ensemble');
%ens_fname = '../work/Prior_Diag.nc'

times = getnc(ens_fname, 'time');
num_times = size(times, 1);

% Select field to plot (u, v, t, ...)
field_name = input('Input field type, 1=u1, 2=v1, 3=t1, 4=u, 5=v, 6=t, 7=qnH, 8=qnOH, 9=qnO')
ens_vec = getnc(ens_fname, field_name);
state_vec = getnc(fname, field_name);

% Ensemble size is
ens_size = size(ens_vec, 2);
                                                                                                               
% Select x and y coordinates
x_coord = input('Select x coordinate');
y_coord = input('Select y coordinate');
% Get level for free atmosphere fields
field_level = input('Input level');

figure(2);
close;
figure(2);
hold on;

% WARNING: MAKE SURE THAT DOING Y THEN X IN RESHAPE IS OKAY
   field = reshape(state_vec(1:num_times,:,:, field_level), [num_times, num_lat, num_lon]);
   plot(times(1:num_times),field(:, y_coord, x_coord), 'r', 'linewidth',2);
   
      ens_mean = reshape(ens_vec(1:num_times,11,:,:, field_level), [num_times num_lat num_lon]);
      ind = ens_mean(:, y_coord, x_coord) < 1.e+36; %% tmp
      plot(times(ind),ens_mean(ind, y_coord, x_coord) ,'k', 'linewidth',2);

% Loop through the ensemble members
   for i = 1 : ens_size -2
      ens = reshape(ens_vec(1:num_times,i,:,:, field_level), [num_times num_lat num_lon]);
      ind = ens(:, y_coord, x_coord) < 1.e+36; %% tmp
      plot(times(ind),ens(ind, y_coord, x_coord));        %% tmp
   end
  
      %ens_spread = reshape(ens_vec(:,12,:,:, field_level), [num_times num_lat num_lon]);
      %ind = ens_spread(:, y_coord, x_coord) < 1.e+36; %% tmp
      %plot(ens_spread(ind, y_coord, x_coord),'g');

title('Time Series of Temperature at lon:90 lat:72.5 height:40km','FontName', 'times','Fontsize',13)
ylabel('Perturbation Temperature [K]','FontName', 'times','Fontsize',13)
xlabel('Time [day]','FontName', 'times','Fontsize',13)
legend('Truth','Ens Mean','Ens Members')
