function state_vec = get_state_copy(fname, copyindex)
%GET_STATE_COPY  Gets a particular copy (one ensemble member) of state from netcdf file
%
% Retrieves a particular copy of a state vector from a file whose
% full or relative path is specified in the file argument.
% NEED TO DEAL WITH ERRORS.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% Need to get a copy with the label copy
copy_meta_data = getnc(fname, 'CopyMetaData');

% Get some information from the truth_file
f = netcdf(fname);
model      = f.model(:);
num_vars   = ncsize(f('StateVariable')); % determine # of state variables
num_copies = ncsize(f('copy')); % determine # of ensemble members
num_times  = ncsize(f('time')); % determine # of output times
close(f);

% Get only the appropriate copy of the state and return
% Should have an error check for bad indices

state_vec = getnc(fname, 'state', [-1, copyindex, -1], [-1, copyindex, -1]);
