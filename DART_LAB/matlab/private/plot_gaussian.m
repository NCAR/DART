function [plot_handle, xlim, ylim] = plot_gaussian(mymean, sd, weight)
%% plot_gaussian Plot gaussian over 5 standard deviations

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

x_min = mymean - 5*sd;
x_max = mymean + 5*sd;
x_range = x_max - x_min;

% Number of points is 1001
num_points = 1001;
interval = x_range / num_points;
x = x_min:interval:x_max;
y = weight * norm_pdf(x, mymean, sd);

plot_handle = plot(x, y);
xlim        = [x_min x_max];
ylim        = [min(y) max(y)];

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
