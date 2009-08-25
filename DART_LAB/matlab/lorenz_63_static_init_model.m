% Initializes class data for L63, sets up global storage
% and reads in control data from input file

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

% Lorenz-63 model parameters
global SIGMA
global R
global B
global DELTAT
global TIME_STEP_DAYS
global TIME_STEP_SECONDS

% Set default values for the model parameters
SIGMA = 10;
R = 28;
B = 8/3;
DELTAT = 0.01;
TIME_STEP_DAYS = 0;
TIMES_STEP_SECONDS = 0;


% Lorenz-63 fixed model parameters
global MODEL_SIZE

MODEL_SIZE = 3;

global STATE_LOC

STATE_LOC = (0:2) / 3;
