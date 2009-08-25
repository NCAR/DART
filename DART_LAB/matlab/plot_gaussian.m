function[plot_handle] = plot_gaussian(mean, sd, weight)
% Plot gaussian over 5 standard deviations

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
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
