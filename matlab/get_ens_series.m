function ens = get_ens_series(fname, state_var_index)
%GET_ENS_SERIES: Returns matrix of time series for all members of ensemble for a variable
%
% the rows of the matrix correspond to time,
% the columns of the matrix correspond to ensemble members
%
% fname = 'Prior_Diag.nc';
% state_var_index = 3;
% ens = get_ens_series(fname,state_var_index);
%

f = netcdf(fname);
model      = f.model(:);
num_vars   = ncsize(f{'StateVariable'}); % determine # of state variables
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
num_times  = ncsize(f{'time'}); % determine # of output times
close(f);

% disp(sprintf('get_ens_series: fname is %s',fname))
% disp(sprintf('get_ens_series: state_var_index is %d',state_var_index))
% disp(sprintf('get_ens_series: num_times   is %d',num_times))
% disp(sprintf('get_ens_series: num_copies  is %d',num_copies))
% disp(sprintf('get_ens_series: num_vars    is %d',num_vars))

metadata    = getnc(fname,'CopyMetaData');           % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's
if ( isempty(copyindices) )
   disp(sprintf('%s has no valid ensemble members',fname))
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   disp(sprintf('%s claims to have %d copies',fname, num_copies))
   error('netcdf file has no ensemble members.')
end
ens_num     = length(copyindices);

% Get the whole thing and then return the ones we want.
% This is usually not too bad, as there are usually many more
% ensemble members than "mean" and "spread" (the two members
% we are NOT interested in for this function).

state_vec = getnc(fname,'state', [-1, -1, state_var_index], ...
                                 [-1, -1, state_var_index]);
% getnc always squeezes out the singleton last dimension.
ens       = state_vec(:,copyindices);

%----------------------------------------------------------------------
% This was way too slow.
%----------------------------------------------------------------------
% Try to loop through each possible ensemble member.
% As long as we keep generating valid copy_indexes, we keep going.
%
%ens_num = 1;
%copy_index = get_copy_index(fname, 'ensemble member1');
%while copy_index > 0,
%   ens(:, ens_num) = get_var_series(fname, copy_index, state_var_index);
%   ens_num = ens_num + 1;    % advance to next potential member
%   copy_string = ['ensemble member', num2str(ens_num)];
%   copy_index = get_copy_index(fname, copy_string);
%end
%ens_num = ens_num - 1;

disp(sprintf('Read %d ensemble members for variable %d in %s', ...
             ens_num, state_var_index,fname));
