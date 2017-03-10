%% state_diag.m uses read_state.m
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

file_name = input('What is file for true state');

true_state = read_state(file_name);

file_name = input('What is file for prior state');

prior_state = read_state(file_name);

file_name = input('What is file for posterior state')

posterior_state = read_state(file_name);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
