function vars = CheckModel(fname);
% CheckModel   tries to ensure that a netcdf file has what we expect. 
%
%

if ( exist(fname) ~= 2 )
   error(sprintf('%s does not exist.',fname))
end

% Get some information from the file
ft         = netcdf(fname);
model      = ft.model(:);
num_vars   = ncsize(ft{'StateVariable'}); % determine # of state variables
num_copies = ncsize(ft{'copy'}); % determine # of ensemble members
num_times  = ncsize(ft{'time'}); % determine # of output times
close(ft); 

if (isempty(model)) 
   error(sprintf('%s has no ''model'' global attribute.',fname))
end
if (prod(size(num_vars)) > 1 ) 
   error(sprintf('%s has no ''StateVariable'' dimension.',fname))
end
if (prod(size(num_copies)) > 1 ) 
   error(sprintf('%s has no ''copy'' dimension.',fname))
end
if (prod(size(num_times)) > 1 ) 
   error(sprintf('%s has no ''time'' dimension.',fname))
end

StateVariable = getnc(fname,'StateVariable');
copy          = getnc(fname,'copy');


switch lower(model)

   case '9var'

      % to be effective for use as a subscript array, 
      % def_state_vars must be a row vector.

      def_state_vars = zeros(1,num_vars);
      def_state_vars(:) = StateVariable(:);

   case 'lorenz_63'

      % to be effective for use as a subscript array, 
      % def_state_vars must be a row vector.

      def_state_vars = zeros(1,num_vars);
      def_state_vars(:) = StateVariable(:);

   case 'lorenz_96'

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_state_vars = round([1 , num_vars/3 , 2*num_vars/3]);

   otherwise

      error(sprintf('model %s unknown',model))

end

vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy), ...
              'def_state_vars',def_state_vars);
