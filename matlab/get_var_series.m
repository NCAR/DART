function var_vec = get_var_series(fname, copynum, state_var)
%GET_VAR_SERIES Gets a particular copy of a state variable from netcdf file
%
% Retrieves a particular copy of a state variable from a file whose
% full or relative path is specified in the file argument.
% NEED TO DEAL WITH ERRORS.
%
% Example 1:
% fname     = '../work/Prior_Diag.nc';
% copynum   = 8;        % Ensemble Member
% state_var = 3;        % which state variable
% var_vec   = get_var_series(fname, copynum, state_var);

f = netcdf(fname);
model      = f.model(:);
num_vars   = ncsize(f{'StateVariable'}); % determine # of state variables
num_copies = ncsize(f{'copy'}); % determine # of ensemble members

% disp(sprintf('get_var_series fname is %s',fname))
% disp(sprintf('get_var_series copynum is %d',copynum))
% disp(sprintf('get_var_series state_var is %d',state_var))
% disp(sprintf('get_var_series model is %s',model))
% disp(sprintf('get_var_series num_vars is %d',num_vars))
% disp(sprintf('get_var_series num_copies is %d',num_copies))

if (copynum > num_copies) 
   disp( sprintf('%s only has %d ''copies/Ensemble members''',fname,num_copies))
   error(sprintf('you wanted copy %d ', copynum))
end

if (state_var > num_vars) 
   disp( sprintf('%s only has %d state variables',fname,num_vars))
   error(sprintf('you wanted variable %d ', state_var))
end

% Get only the appropriate copy of the state and return
var_vec = getnc(fname, 'state', [-1, copynum, state_var], ...
                                [-1, copynum, state_var]);

