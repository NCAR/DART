%% state_diag.m uses read_state.m
%

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

file_name = input('What is file for true state');

true_state = read_state(file_name);

file_name = input('What is file for prior state');

prior_state = read_state(file_name);

file_name = input('What is file for posterior state')

posterior_state = read_state(file_name);

