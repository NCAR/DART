function list_copies(fname)
% LIST COPIES  Lists the copy_meta_data for an input relative or absolute path.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 

fprintf('------------------------------------------\n');
fprintf('Contents of netcdf file %s\n ', fname);

%file_title = getnc(fname);
%fprintf('Title is %s ', file_title);

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
   fprintf('Copy %3i is %s \n', i, copy_meta_data(i, :)) ;
end
fprintf('------------------------------------------\n');



