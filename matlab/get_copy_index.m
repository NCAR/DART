function copy_index = get_copy_index(fname, copystring)
%GET_COPY_INDEX  Gets an index corresponding to copy meta_data string
%
% Retrieves index associated with a given copy meta_data string in 
% file fname. If string is not found in meta_data list, a -1 is returned.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

f = netcdf(fname);
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
close(f)

% Need to get the Copy MetaData (strings with the names)
% We then search the metadata for each copy for the appropriate copies.
copy_meta_data = getnc(fname, 'CopyMetaData');

% For a single copy, the size is nx1, for more k copies, it's kxn
if size(copy_meta_data, 2) == 1
   copy_meta_data = transpose(copy_meta_data);
end

% Figure out which copy is the matching one
copy_index = -1;
for i = 1:num_copies,
   if strcmp(deblank(copy_meta_data(i, :)) , copystring) == 1
%      fprintf('copy %3i selected \n', i);
      copy_index = i;
   end
end

% Provide modest error support

if (copy_index < 0)
   disp(sprintf('ERROR: %s is not a valid metadata string for file %s', ...
                copystring, fname))
   disp('valid metadata strings are: ')
   for i = 1:num_copies,
      disp(sprintf('%s',deblank(copy_meta_data(i,:))))
   end
   error('please try again')
end
