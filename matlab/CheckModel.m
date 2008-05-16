function vars = CheckModel(fname);
% CheckModel   tries to ensure that a netcdf file has what we expect. 
%
% vars is a structure containing a minimal amount of metadata about the netCDF file.
% 
% EXAMPLE:
% fname = 'Prior_Diag.nc';
% vars = CheckModel(fname) 

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
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

   case {'9var','lorenz_63','lorenz_84','ikeda'}

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

   case 'forced_lorenz_96'

      % This model has the state variables replicated, so there is a difference
      % between num_state_vars and the length of the state variable.
      forcing           = f.model_forcing(:);
      delta_t           = f.model_delta_t(:);
      time_step_days    = f.model_time_step_days(:);
      time_step_seconds = f.model_time_step_seconds(:);
      num_model_vars    = cast(f.model_num_state_vars(:),'double');  % ACTUAL state vars, not

      atts = dim(f{'StateVariable'}); num_vars  = length(atts{1}); % # of state varbls
      if (prod(size(num_vars)) > 1 ) 
         error(sprintf('%s has no ''StateVariable'' dimension.',fname))
      end
      StateVariable = getnc(fname,'StateVariable');

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_state_vars = round([1 , num_model_vars/3 , 2*num_model_vars/3]);
      def_force_vars = num_model_vars + def_state_vars;

      vars = struct('model',model, ...
              'def_var','state', ...
              'num_state_vars',num_vars, ...
              'num_model_vars',num_model_vars, ...
              'num_force_vars',num_vars - num_model_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',      1, 'max_state_var', num_vars, ...
              'min_model_var',      1, 'max_model_var', num_model_vars, ...
              'min_force_var',      1, 'max_force_var', num_vars - num_model_vars, ...
              'min_ens_mem',min(copy), 'max_ens_mem',   max(copy), ...
              'def_state_vars',def_state_vars, ...
              'def_force_vars',def_force_vars);

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

   case 'simple_advection'

      loc1d = getnc(fname,'loc1d');
      num_locs = length(loc1d);

      if ( isempty(f{'state'}) )
         varnames = {'concentration','source','wind', ...
                    'mean_source','source_phase'};
         def_inds = round([1 , num_locs/3 , 2*num_locs/3]);
      else
         varnames = {'state'};
         def_inds = [1 13 27];
      end

      vars = struct('model'       ,model, ...
              'loc1d'             ,loc1d, ...
              'num_ens_members'   ,num_copies, ...
              'min_ens_mem'       ,min(copy), ...
              'max_ens_mem'       ,max(copy), ...
              'time_series_length',num_times, ...
              'model_size'        ,length(varnames)*length(loc1d), ...
              'def_var'           ,varnames{1}, ...
              'min_state_var'     ,1, ...
              'max_state_var'     ,num_locs, ...
              'def_state_vars'    ,def_inds, ...
              'num_vars'          ,length(varnames));
      vars.vars = varnames;

   case 'fms_bgrid'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);      % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)

      varnames  = {'ps','t','u','v'};
      num_vars = length(varnames);
      nlevels   = ncsize(f('lev')); % determine # of state variables
      if (prod(size(nlevels)) > 1 ) 
         error(sprintf('%s has no ''lev'' dimension.',fname))
      end
%     times  = getnc(fname,'time');
%     TmpI   = getnc(fname,'TmpI');    % longitude
%     TmpJ   = getnc(fname,'TmpJ');    % latitude
%     levels = getnc(fname,'level');
%     VelI   = getnc(fname,'VelI');    % longitude
%     VelJ   = getnc(fname,'VelJ');    % latitude

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy));

      vars.vars = varnames;

   case 'cam'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);     % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)
      % have not yet figured out a way to only use non-coordinate variables.

      varnames = {'PS','T','U','V','Q','CLDLIQ','CLDICE'};
      num_vars = length(varnames);
      nlevels  = ncsize(f('lev')); % determine # of state variables
      if (prod(size(nlevels)) > 1 ) 
          error(sprintf('%s has no ''lev'' dimension.',fname))
      end

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );
         %    'max_ens_mem',max(copy), ...
         %    'varnames',varnames);

   case 'pbl_1d'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);     % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)

      num_vars  = 22; % ps, t, u, v
      z_level   = ncsize(f('z_level')); % determine # of state variables
      sl_level  = ncsize(f('sl_level')); % determine # of state variables
      if (prod(size(z_level)) > 1 ) 
         error(sprintf('%s has no ''z_level'' dimension.',fname))
      end
      times    = getnc(fname,'time');
      z_level  = getnc(fname,'z_level');
      sl_level = getnc(fname,'sl_level');

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy));

   case 'pe2lyr'

      % Since this is a 3D model, only the most rudimentary information
      % is gathered here. Each plot requires different information,
      % so there is a separate function (GetPe2lyrInfo.m) that gets
      % the information for each specific plot type.

      varnames  = {'u','v','z'};
      num_vars  = length(varnames);

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );

      vars.vars = varnames;

   case 'mitgcm_ocean'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);     % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)
      % have not yet figured out a way to only use non-coordinate variables.

      varnames = {'S','T','U','V','SSH'};
      num_vars = length(varnames);
      nlevels  = ncsize(f('ZG')); % determine # of state variables
      if (prod(size(nlevels)) > 1 ) 
          error(sprintf('%s has no ''ZG'' dimension.',fname))
      end

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );

       vars.vars = varnames;

   otherwise

      error(sprintf('model %s unknown',model))

end
close(f); 

