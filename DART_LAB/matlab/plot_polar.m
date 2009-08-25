function h = plot_polar(y, x, mean_dist, string, model_size)
%

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

% Y includes a wraparound point, x does not
x_t(model_size + 1) = x(1);
x_t(1:model_size) = x;
h = polar(y, mean_dist + x_t, string);

end
