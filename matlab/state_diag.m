% state_diag.m uses read_state.m
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2006, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
file_name = input('What is file for true state');

true_state = read_state(file_name);

file_name = input('What is file for prior state');

prior_state = read_state(file_name);

file_name = input('What is file for posterior state')

posterior_state = read_state(file_name);

