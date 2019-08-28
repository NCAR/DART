function L63 = lorenz_63_static_init_model()

%% lorenz_63_static_init_model Initializes class data for L63, sets up global storage
% and reads in control data from input file

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Lorenz-63 model parameters

L63.sigma             = 10;
L63.r                 = 28;
L63.b                 = 8/3;
L63.deltat            = 0.01;
L63.time_step_days    = 0;
L63.time_step_seconds = 0;
L63.model_size        = 3;
L63.state_loc         = (0:2) / 3;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
