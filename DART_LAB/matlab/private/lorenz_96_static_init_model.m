%% lorenz_96_static_init_model Initializes class data for L96, sets up global storage
% and reads in control data from input file

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Lorenz-96 model parameters
global FORCING
global DELTA_T
global TIME_STEP_DAYS
global TIME_STEP_SECONDS

% Set default values for the model parameters
FORCING = 8;
DELTA_T = 0.05;
TIME_STEP_DAYS = 0;
TIME_STEP_SECONDS = 0;


% Lorenz-96 fixed model parameters
global MODEL_SIZE

MODEL_SIZE = 40;

global STATE_LOC

STATE_LOC = (0:MODEL_SIZE - 1) / MODEL_SIZE;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

