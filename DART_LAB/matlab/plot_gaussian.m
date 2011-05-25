function[plot_handle] = plot_gaussian(mean, sd, weight)
%% plot_gaussian Plot gaussian over 5 standard deviations

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

x_min = mean - 5*sd;
x_max = mean + 5*sd;
x_range = x_max - x_min;

% Number of points is 1001
num_points = 1001;
interval = x_range / num_points;
x = x_min:interval:x_max;
y = weight * normpdf(x, mean, sd);

plot_handle = plot(x, y);
