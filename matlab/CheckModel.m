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
f         = netcdf(fname);
model     = f.model(:);

atts = dim(f{'copy'}); num_copies = length(atts{1}); % determine # of ensemble members
atts = dim(f{'time'}); num_times  = length(atts{1}); % determine # of output times

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

   case {'9var','lorenz_63','lorenz_84'}

      atts = dim(f{'StateVariable'}); num_vars  = length(atts{1}); % # of state varbls
      if (prod(size(num_vars)) > 1 ) 
         error(sprintf('%s has no ''StateVariable'' dimension.',fname))
      end
      StateVariable = getnc(fname,'StateVariable');

      def_state_vars = zeros(1,num_vars);    % for use as a subscript array, 
      def_state_vars(:) = StateVariable(:);  % def_state_vars must be a row vector.

      vars = struct('model',model, ...
              'def_var','state', ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy), ...
              'def_state_vars',def_state_vars);

   case {'lorenz_96', 'lorenz_04'}

      atts = dim(f{'StateVariable'}); num_vars  = length(atts{1}); % # of state varbls
      if (prod(size(num_vars)) > 1 ) 
         error(sprintf('%s has no ''StateVariable'' dimension.',fname))
      end
      StateVariable = getnc(fname,'StateVariable');

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_state_vars = round([1 , num_vars/3 , 2*num_vars/3]);

      vars = struct('model',model, ...
              'def_var','state', ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy), ...
              'def_state_vars',def_state_vars);

   case 'lorenz_96_2scale'

      atts = dim(f{'Xdim'}); num_X  = length(atts{1}); % # of X variables
      if (prod(size(num_X)) > 1 ) 
         error(sprintf('%s has no ''Xdim'' dimension.',fname))
      end
      Xdim = getnc(fname,'Xdim');

      atts = dim(f{'Ydim'}); num_Y  = length(atts{1}); % # of Y variables
      if (prod(size(num_Y)) > 1 ) 
         error(sprintf('%s has no ''Ydim'' dimension.',fname))
      end
      Ydim = getnc(fname,'Ydim');

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_X_inds = round([1 , num_X/3 , 2*num_X/3]);
      def_Y_inds = round([1 , num_Y/3 , 2*num_Y/3]);

      vars = struct('model',model, ...
              'def_var','X', ...
              'num_state_vars',num_X, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',min(Xdim), 'max_state_var',max(Xdim), ...
              'min_X_var',    min(Xdim), 'max_X_var',    max(Xdim), ...
              'min_Y_var',    min(Ydim), 'max_Y_var',    max(Ydim), ...
              'min_ens_mem',  min(copy), 'max_ens_mem',  max(copy), ...
              'def_state_vars',def_X_inds);

   case 'fms_bgrid'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);     % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)

      num_vars  = 4; % ps, t, u, v
      nlevels   = ncsize(f('lev')); % determine # of state variables
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

   case 'pe2lyr'

   otherwise

      error(sprintf('model %s unknown',model))

end
close(f); 

