%GET_COPY_INDEX  Gets a the index corresponding to copy meta_data string

% Retrieves index associated with a given copy meta_data string in 
% file fname. If string is not found in meta_data list, a -1 is returned.

function index = get_copy_index(copy, fname)

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
index = -1;
for i = 1:num_copies
   if strcmp(deblank(copy_meta_data(i, :)) , copy) == 1
%      fprintf('copy %3i selected \n', i);
      index = i;
   end
end
