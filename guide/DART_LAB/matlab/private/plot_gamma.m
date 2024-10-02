function [plot_handle] = plot_gamma(shape, scale)
%% plot_gamma Plot gamma over 5 standard deviations

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$


% Get the mean and variance for this shape and scale
[gmean, gvar] = gamstat(shape, scale);
gsd = sqrt(gvar);

% Get points to a range where pdf is small
x_min = 0;
x_max = gmean + 5*gsd;

% Number of points is 1001
num_points = 1001;
interval = x_max / num_points;
x = 0:interval:x_max;
y = gampdf(x, shape, scale);

plot_handle = plot(x, y);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
