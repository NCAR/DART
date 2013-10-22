function h = plot_polar(y, x, mean_dist, string, model_size)
%% plot_polar

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Y includes a wraparound point, x does not
x_t(model_size + 1) = x(1);
x_t(1:model_size) = x;
h = polar(y, mean_dist + x_t, string);

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

