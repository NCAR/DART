function ens = get_ens_series(fname, state_var_index)
%GET_ENS_SERIES: Returns array of time series for all members of ensemble for varible

f = netcdf(fname);
model      = f.model(:);
num_vars   = ncsize(f{'StateVariable'}); % determine # of state variables
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
num_times  = ncsize(f{'time'}); % determine # of output times

% disp(sprintf('get_ens_series: fname is %s',fname))
% disp(sprintf('get_ens_series: state_var_index is %d',state_var_index))
% disp(sprintf('get_ens_series: num_times   is %d',num_times))
% disp(sprintf('get_ens_series: num_copies  is %d',num_copies))
% disp(sprintf('get_ens_series: num_vars    is %d',num_vars))

% Try to loop through each possible ensemble member.
% As long as we keep generating valid copy_indexes, we keep going.

% It would be more efficient to search though the metadata 
% to collect the useful copy indices and then preallocate
% storage ...

ens_num = 1;
copy_index = get_copy_index(fname, 'ensemble member1');
while copy_index > 0,
   ens(:, ens_num) = get_var_series(fname, copy_index, state_var_index);
   ens_num = ens_num + 1;    % advance to next potential member
   copy_string = ['ensemble member', num2str(ens_num)];
   copy_index = get_copy_index(fname, copy_string);
end

% disp(sprintf('Total number of ensemble members is %d', ens_num - 1));
disp(sprintf('Read %d ensemble members for state variable %d', ...
             ens_num-1, state_var_index));
