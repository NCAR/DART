function vars = CheckModel(fname);
% CheckModel   tries to ensure that a netcdf file has what we expect. 

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% $Source$
% $Revision$
% $Date$

if ( exist(fname) ~= 2 )
   error(sprintf('%s does not exist.',fname))
end

% Get some information from the file
ft         = netcdf(fname);
model      = ft.model(:);
num_copies = ncsize(ft('copy')); % determine # of ensemble members
num_times  = ncsize(ft('time')); % determine # of output times

if (isempty(model)) 
   error(sprintf('%s has no ''model'' global attribute.',fname))
end
if (prod(size(num_copies)) > 1 ) 
   error(sprintf('%s has no ''copy'' dimension.',fname))
end
if (prod(size(num_times)) > 1 ) 
   error(sprintf('%s has no ''time'' dimension.',fname))
end

copy          = getnc(fname,'copy');

switch lower(model)

   case '9var'

      num_vars   = ncsize(ft{'StateVariable'}); % determine # of state variables
      if (prod(size(num_vars)) > 1 ) 
         error(sprintf('%s has no ''StateVariable'' dimension.',fname))
      end
      StateVariable = getnc(fname,'StateVariable');

      def_state_vars = zeros(1,num_vars);    % for use as a subscript array, 
      def_state_vars(:) = StateVariable(:);  % def_state_vars must be a row vector.

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy), ...
              'def_state_vars',def_state_vars);

   case 'lorenz_63'

      num_vars   = ncsize(ft{'StateVariable'}); % determine # of state variables
      if (prod(size(num_vars)) > 1 ) 
         error(sprintf('%s has no ''StateVariable'' dimension.',fname))
      end
      StateVariable = getnc(fname,'StateVariable');

      def_state_vars = zeros(1,num_vars);    % for use as a subscript array, 
      def_state_vars(:) = StateVariable(:);  % def_state_vars must be a row vector.

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy), ...
              'def_state_vars',def_state_vars);

   case 'lorenz_96'

      num_vars   = ncsize(ft{'StateVariable'}); % determine # of state variables
      if (prod(size(num_vars)) > 1 ) 
         error(sprintf('%s has no ''StateVariable'' dimension.',fname))
      end
      StateVariable = getnc(fname,'StateVariable');

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_state_vars = round([1 , num_vars/3 , 2*num_vars/3]);

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy), ...
              'def_state_vars',def_state_vars);

   case 'fms_bgrid'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(ft);     % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)

      num_vars  = 4; % ps, t, u, v
      nlevels   = ncsize(ft('lev')); % determine # of state variables
      if (prod(size(nlevels)) > 1 ) 
         error(sprintf('%s has no ''lev'' dimension.',fname))
      end
      times  = getnc(fname,'time');
      TmpI   = getnc(fname,'TmpI');    % longitude
      TmpJ   = getnc(fname,'TmpJ');    % latitude
      levels = getnc(fname,'level');
      VelI   = getnc(fname,'VelI');    % longitude
      VelJ   = getnc(fname,'VelJ');    % latitude

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy));

   otherwise

      error(sprintf('model %s unknown',model))

end
close(ft); 

