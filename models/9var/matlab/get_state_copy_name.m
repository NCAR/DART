function state_vec = get_state_copy_name(copy, fname)
% GET_STATE_COPY_NAME  Gets a particular copy of state from netcdf file
%
% Retrieves a particular copy of a state vector from a file whose
% full or relative path is specified in the file argument.
% NEED TO DEAL WITH ERRORS.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
loc = getnc(fname, 'loc1d');
num_loc = size(loc, 1);
times = getnc(fname, 'time');
num_times = size(times, 1);

% Need to get a copy with the label copy
copy_meta_data = getnc(fname, 'CopyMetaData');
% For a single copy, the size is nx1, for more k copies, it's kxn
if size(copy_meta_data, 2) == 1
   copy_meta_data = transpose(copy_meta_data);
end
num_copies = size(copy_meta_data, 1);

% Figure out which copy is the matching one
index = 0;
for i = 1:num_copies
   fprintf('meta data for copy %3i is %s \n', i, copy_meta_data(i, :)) ;
   if strncmp(copy_meta_data(i, :) , copy, size(copy, 2)) == 1
      fprintf('copy %3i selected \n', i);
      index = i;
   end
end

if index ==0 
   fprintf('ERROR: copy %s not found in file \n', copy);
   state_vec = [];
   return
end

% Get only the appropriate copy of the state and return
if num_copies == 1
   state_vec = getnc(fname, 'state');
else
   state_vec = getnc(fname, 'state', [-1, index, -1], [-1, index, -1]);
end

