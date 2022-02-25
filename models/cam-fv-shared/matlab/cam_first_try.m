%% cam_first_try

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Get file name of true state fileS
fname = input('Input true state name');
%fname = 'True_State.nc';
lon = getnc(fname, 'I');
num_lon = size(lon, 1);
lat = getnc(fname, 'J');
num_lat = size(lat, 1);
level = getnc(fname, 'level');
num_level = size(level, 1);

state_vec = getnc(fname, 'state');

% Get a time level from the user
time_ind = input('Input time level');

single_state = state_vec(time_ind, :);

% Select field to plot (ps, t, u, v)
field_num = input('Input field type, 1=ps, 2=t, 3=u, or 4=v, 5=q')

% Get level for free atmosphere fields
if field_num > 1
   field_level = input('Input level');
else
   field_level = 1;
end

% Extract one of the fields
   offset = (field_level - 1) * 4 + field_num
% Stride is number of fields in column
   stride = num_level * 4 + 1
   field_vec = single_state(offset : stride : (stride) * (num_lon * num_lat));
   field = reshape(field_vec, [num_lat, num_lon]); 


[C, h] = contourf(field);
clabel(C, h);

% Loop for another try
cam_first_try;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
