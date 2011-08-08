function vars = CheckModel(fname);
%% CheckModel   tries to ensure that a netcdf file has what we expect. 
%
% vars is a structure containing a minimal amount of metadata about the netCDF file.
% 
% EXAMPLE:
% fname = 'Prior_Diag.nc';
% vars = CheckModel(fname) 

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

% Get some information from the file
model  = nc_attget(fname,nc_global,'model');

num_copies = dim_length(fname,'copy'); % determine # of ensemble members
num_times  = dim_length(fname,'time'); % determine # of output times

if (isempty(model)) 
   error('%s has no ''model'' global attribute.',fname)
end

copy = nc_varget(fname,'copy');

switch lower(model)

   case {'9var','lorenz_63','lorenz_84','ikeda'}

      num_vars      = dim_length(fname,'StateVariable'); % determine # of state varbls
      StateVariable =  nc_varget(fname,'StateVariable');

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

      vars.fname = fname;

   case {'lorenz_96', 'lorenz_04'}

      num_vars      = dim_length(fname,'StateVariable'); % determine # of state varbls
      StateVariable =  nc_varget(fname,'StateVariable');

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

      vars.fname = fname;

   case 'forced_lorenz_96'

      % This model has the state variables replicated, so there is a difference
      % between num_state_vars and the length of the state variable.
      forcing           = nc_attget(fname, nc_global, 'model_forcing');
      delta_t           = nc_attget(fname, nc_global, 'model_delta_t');
      time_step_days    = nc_attget(fname, nc_global, 'model_time_step_days');
      time_step_seconds = nc_attget(fname, nc_global, 'model_time_step_seconds');
      num_model_vars    = nc_attget(fname, nc_global, 'model_num_state_vars');

      num_vars      = dim_length(fname,'StateVariable'); % determine # of state varbls
      StateVariable =  nc_varget(fname,'StateVariable');

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

      vars.fname = fname;

   case 'lorenz_96_2scale'

      num_X  = dim_length(fname,'Xdim'); % # of X variables
      Xdim   =  nc_varget(fname,'Xdim');

      num_Y  = dim_length(fname,'Ydim'); % # of Y variables
      Ydim   =  nc_varget(fname,'Ydim');

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

      vars.fname = fname;

   case 'simple_advection'

      num_locs = dim_length(fname,'loc1d'); % # of X variables
      loc1d    =  nc_varget(fname,'loc1d');

      if ( nc_isvar(fname,'state') )
         varnames = {'state'};
         def_inds = [1 13 27];
      else
         varnames = {'concentration','source','wind', ...
                    'mean_source','source_phase'};
         def_inds = round([1 , num_locs/3 , 2*num_locs/3]);
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

      vars.vars  = varnames;
      vars.fname = fname;

   case 'fms_bgrid'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);      % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)

      varnames = {'ps','t','u','v'};
      num_vars = length(varnames);
      nlevels  = dim_length(fname,'lev'); % determine # of state variables

%     times  = nc_varget(fname,'time');
%     TmpI   = nc_varget(fname,'TmpI');    % longitude
%     TmpJ   = nc_varget(fname,'TmpJ');    % latitude
%     levels = nc_varget(fname,'level');
%     VelI   = nc_varget(fname,'VelI');    % longitude
%     VelJ   = nc_varget(fname,'VelJ');    % latitude

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy));

      vars.vars  = varnames;
      vars.fname = fname;

   case {'cam','tiegcm'}

      varnames = get_DARTvars(fname);
      num_vars = length(varnames);
      nlevels  = dim_length(fname,'lev'); % determine # of state variables

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );

      vars.vars  = varnames;
      vars.fname = fname;

   case 'pbl_1d'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);     % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)

      num_vars = 22; % ps, t, u, v
      z_level  = dim_length(fname, 'z_level'); % determine # of state variables
      sl_level = dim_length(fname,'sl_level'); % determine # of state variables
      times    = nc_varget(fname,'time');
      z_level  = nc_varget(fname,'z_level');
      sl_level = nc_varget(fname,'sl_level');

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy));

      vars.fname = fname;

   case 'pe2lyr'

      % Since this is a 3D model, only the most rudimentary information
      % is gathered here. Each plot requires different information,
      % so there is a separate function (GetPe2lyrInfo.m) that gets
      % the information for each specific plot type.

      varnames = {'u','v','z'};
      num_vars = length(varnames);

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );

      vars.vars  = varnames;
      vars.fname = fname;

   case 'mitgcm_ocean'

      % A more robust way would be to use the netcdf low-level ops:
      % bob = var(f);     % bob is a cell array of ncvars
      % name(bob{1})       % is the variable name string
      % bob{1}(:)          % is the value of the netcdf variable  (no offset/scale)
      % have not yet figured out a way to only use non-coordinate variables.

      varnames = {'S','T','U','V','SSH'};
      num_vars = length(varnames);
      nlevels  = dim_length(fname,'ZG'); % determine # of state variables

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );

      vars.vars  = varnames;
      vars.fname = fname;

   case 'wrf'

      % requires a 'domain' and 'bottom_top_d01' dimension.
      % without both of these, it will fail in an ugly fashion.

      varnames    = get_DARTvars(fname);
      num_vars    = length(varnames);
      dinfo       = nc_getdiminfo(fname,'domain');
      num_domains = dinfo.Length;
      dinfo       = nc_getdiminfo(fname,'bottom_top_d01');
      num_levels  = dinfo.Length;

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'num_unstaggered_levels',num_levels, ...
              'num_domains',num_domains, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy));

      vars.vars  = varnames;
      vars.fname = fname;

   otherwise

      error('model %s unknown',model)

end




function x = dim_length(fname,dimname)

y = nc_isvar(fname,dimname);
if (y < 1)
   error('%s has no %s dimension/coordinate variable',fname,dimname)
end
bob = nc_getdiminfo(fname,dimname);
x   = bob.Length;

