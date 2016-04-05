% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
file_name = input('What is file for true state');

read_state;

true_state = state;

file_name = input('What is file for prior state');

read_state;

prior_state = state;

file_name = input('What is file for posterior state')

read_state;

posterior_state = state;

