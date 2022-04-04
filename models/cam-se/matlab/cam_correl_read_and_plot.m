%% cam_correl_read_and_plot
% usage: assumes cam_correl has already been run and has set a slew of global
% variables.  these can be changed and then this subroutine rerun again.  
% if you are changing the input files you must rerun this entire routine.  
% if you are just changing plotting parameters, you can rerun 'cam_correl_plot'
% only.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% doc for Get_Field, Get_Point funcs syntax
% x[nlat,nlon,nens] = Get_Field(filename,basevar,timeind,level,ens_size);  
% x[nens] = Get_Point(filename,basevar,timeind,lat,lon,level,ens_size);  
%
% some quick sanity checks first, and then some setup before looping
% anything which changes based on the input files changing must be here.
% anything which relates only to the plot or probe points is in the 
% plot-only routine.
%


% assume all files are conformant, so pick one and use the parms from it.
filename = sprintf(file_base_pattern, file_ids(1));

% how many timesteps per file?
tt = getnc(filename,'time');
steps_per_file = length(tt);

% grid resolution of data
lat = getnc(filename,'lat');
lat_res = length(lat);
lon = getnc(filename,'lon');
lon_res = length(lon);
lev = getnc(filename,'lev');
lev_res = length(lev);

% should match number of file ids specified above
nfiles = floor((num_lags + steps_per_file-1) / steps_per_file);   
if (nfiles ~= length(file_ids)) 
  disp(sprintf('number of file ids: %d', length(file_ids)));
  disp(sprintf('number of timesteps/file: %d', steps_per_file));
  disp(sprintf('number of lags: %d', num_lags));
  error('number of file_ids must match num_lags * number of timesteps/file')
end
  
% ensemble size: get from num_copies - 2 for data mean/var - 2 for inf mean/var
ncopies = getnc(filename,'copy');
ens_size = length(ncopies) - 4;

% get coordinate variables and metadata from the netCDF file.

% metadata (units, time, etc)
f = netcdf(filename,'nowrite');
base_varunits  = f{base_var}.units(:);
probe_varunits = f{probe_var}.units(:);
timebase = f{'time'}.units(:);
model    = f.model(:);
close(f);


% Make space for x, t
x(1:num_lags, 1:ens_size, 1:lat_res, 1:lon_res) = 0.0;
t(num_lags) = 0;

% read in data and time stamps for all files, all times
lag = 1;
for i = 1:nfiles
  filename = sprintf(file_base_pattern, file_ids(i));
  tt = getnc(filename,'time');

  for timeind = 1:steps_per_file

    x(lag, :, :, :) = Get_Field(filename, base_var, timeind, ...
                                base_level, ens_size);

    t(lag) = tt(timeind);

    lag = lag + 1;

  end
end

%
%  loop and construct the plots
%

cam_correl_plot

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
