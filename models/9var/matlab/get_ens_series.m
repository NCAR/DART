%GET_ENS_SERIES: Returns array of time series for all members of ensemble for varible

function ens = get_ens_series(var, fname)

% Need to get a copy with the label copy
copy_meta_data = getnc(fname, 'CopyMetaData');
% For a single copy, the size is nx1, for more k copies, it's kxn
if size(copy_meta_data, 2) == 1
   copy_meta_data = transpose(copy_meta_data);
end
num_copies = size(copy_meta_data, 1);

% Try to loop through each possible ensemble member
ens_num = 1;
index = get_copy_index('ensemble member1', fname);
while index > 0,
   ens(:, ens_num) = get_var_series(var, index, fname);
   ens_num = ens_num + 1;
   copy_string = ['ensemble member', num2str(ens_num)];
   index = get_copy_index(copy_string, fname);
end

fprintf('Total number of ensemble members is %d \n', ens_num - 1);
