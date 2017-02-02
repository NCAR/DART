%% cam_correl
% usage: computes forward/backwards correlations based on data
% at a single lat,lon.   the script has many global variables to set
% and first calls a 'read data' function (which may be slow depending
% on the size of the diagnostic files) and then a 'plot' function.
% the intent is that after setting the globals and running this
% script once, you can reset the globals to different values and replot
% without waiting to reread the data (assuming the data files and which
% variable being plotted doesn't change.)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% global parms - set once here.   none of these have semicolons so they
%                echo their values as they print.
%

% direction - +1 for forward or -1 for backwards/adjoint
direction = +1

% which variable to use against probe point, and which model level
base_var   = 'T'
base_level = 19

% which lat,lon point and variable field to do the correlations with
probe_var    = 'T'
probe_level  = 19
probe_pt_lat = 100
probe_pt_lon = 160

% limits of the plot area
% x axis:
plot_lon_min = 195
plot_lon_max = 285
plot_lon_ntics = 3    % including the 2 end ticks

% y axis:
plot_lat_min =  30
plot_lat_max =  70
plot_lat_ntics = 5    % including the 2 end ticks

% how much of center of colorbar to white out to get just the extremes
% the full colorbar is 64 bins, centered on 0.  this amount will be
% blanked on either side, so e.g. 2 will give 4 blank bins.
colorbar_blank = 4

% input data files.
% filename -- common patterns:  01_NN/Prior_Diag.nc, Prior_Diag_%d.nc
% use the 'sprintf' syntax - literal chars plus %d for any (integer) numbers.
% put time_steps in increasing/forward order for both forward & adjoint.
file_base_pattern = 'Prior_Diag_%d.nc'
file_ids = [52,53]

% number of time lags forward or backward to cross-correlate
num_lags = 4    % this must match the number of files and ntimesteps/file

% print output?  0=no, 1=yesi.  and if so, what format? 
print_to_file = 0

%print_format = '-deps'     % b/w postscript
%print_format = '-depsc'    % color postscript
 print_format = '-depsc2'   % color postscript, level 2
%print_format = '-depsc2 -tiff' % color ps, with tiff preview
%print_format = '-dpdf'     % color pdf
%print_format = '-dpng'     % color portable network graphics

% for some uses you want the plot with no colorbar and no labels on 
% any of the axis.  for that, set 'nolabels' to 1 (0 prints normal labels)
% this will still print a title since it contains the timestamp which
% varies from plot to plot - but if you don't want the title, set
% 'notitle' to 1.  (0 prints like normal)
nolabels = 0     % no axis labels or colorbars on the plots, title remains
notitle  = 0     % no title 

%
% end of global variable section
%

%
% the file reading and plotting are in subroutines to make it easy to reset the
% global vars and retry again.
%

cam_correl_read_and_plot

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
