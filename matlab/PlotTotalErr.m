function PlotTotalErr(truth_file,diagn_file)
% Plots summary plots of global error and spread
%
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Posterior_Diag.nc';
% PlotTotalErr(truth_file,diagn_file);

if ( exist(truth_file) ~= 2 ), error(sprintf('(truth_file) %s does not exist.',truth_file)), end
if ( exist(diagn_file) ~= 2 ), error(sprintf('(diagn_file) %s does not exist.',diagn_file)), end

% Get some information from the truth_file 
ft = netcdf(truth_file);
t.model      = ft.model(:);
t.num_vars   = ncsize(ft{'StateVariable'}); % determine # of state variables
t.num_copies = ncsize(ft{'copy'}); % determine # of ensemble members
t.num_times  = ncsize(ft{'time'}); % determine # of output times
close(ft);

% Get some information from the diagn_file 
fd = netcdf(diagn_file);
d.model      = fd.model(:);
d.num_vars   = ncsize(fd{'StateVariable'}); % determine # of state variables
d.num_copies = ncsize(fd{'copy'}); % determine # of ensemble members
d.num_times  = ncsize(fd{'time'}); % determine # of output times
close(fd);

% rudimentary bulletproofing
if (strcmp(t.model,d.model) ~= 1)
   error('no No NO ... models must be the same')
end

% Get the netcdf variable indices for desired "copies"
% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(truth_file, 'true state' );
ens_mean_index   = get_copy_index(diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(diagn_file, 'ensemble spread');

% Get some useful plotting arrays
times = getnc(truth_file,'time');

switch lower(t.model)

   case '9var'

      % Get the appropriate netcdf variables
      truth  = get_state_copy(truth_file,     truth_index);
      ens    = get_state_copy(diagn_file,  ens_mean_index);
      spread = get_state_copy(diagn_file,ens_spread_index);

      % Also need to compute the spread; zero truth for this and compute
      % distance from 0 (verify algorithm)
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);

      clf;
      plot(times,err, 'b', times,err_spread, 'r');
      legend('Ensemble Mean Total Error','Ensemble Spread',0)
      title(sprintf('model %s Total Error over all %d variables',...
                     t.model, d.num_vars),'interpreter','none')
      xlabel(sprintf('model time (%d timesteps)',t.num_times))
      ylabel('Error')

   case 'lorenz_63'

      % Get the appropriate netcdf variables
      truth  = get_state_copy(truth_file,     truth_index);
      ens    = get_state_copy(diagn_file,  ens_mean_index);
      spread = get_state_copy(diagn_file,ens_spread_index);

      % Also need to compute the spread; zero truth for this and compute
      % distance from 0 (verify algorithm)
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);

      clf;
      plot(times,err, 'b', times,err_spread, 'r');
      legend('Ensemble Mean Total Error','Ensemble Spread',0)
      title(sprintf('model %s Total Error over all %d variables',...
                     t.model, d.num_vars),'interpreter','none')
      xlabel(sprintf('model time (%d timesteps)',t.num_times))
      ylabel('Total Error')

   case 'lorenz_96'
      disp('lorenz_96 not implemented yet.')

   otherwise
end



