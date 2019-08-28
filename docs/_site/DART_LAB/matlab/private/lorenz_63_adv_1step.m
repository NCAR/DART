function[x_new, time_new] = lorenz_63_adv_1step(x, time)
%% lorenz_63)adv_1step advances the lorenz convective 3 variable model
% for a single two step runge-kutta time step
%
% x is the 3-vector state, time is the 2-vector days and seconds time

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

global DELTAT

% Set the fraction for the rk-2 time step
fract = 1;

% Compute the first intermediate step
dx = comp_dt(x);
x1 = x + fract * DELTAT * dx;

% Compute the second intermediate step
dx = comp_dt(x1);
x2 = x1 + fract * DELTAT * dx;

% New value for x is average of original value and second intermediate
x_new = (x + x2) / 2;

% Update the time ; Non-dimensional single unit for development
time_new = time + 1;

end


function[dt] = comp_dt(x)

global SIGMA
global R
global B

dt(1) = SIGMA * (x(2) - x(1));
dt(2) = -1 * x(1)*x(3) + R*x(1) - x(2);
dt(3) = x(1)*x(2) - B*x(3);

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
