function PlotEnsMeanTimeSeries(truth_file,diagn_file)
% Plots time series of ensemble members, mean and truth
%
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Prior_Diag.nc';
% PlotEnsMeanTimeSeries(truth_file,diagn_file)

if ( exist(truth_file) ~= 2 ), error(sprintf('(truth_file) %s does not exist.',truth_file)), end
if ( exist(diagn_file) ~= 2 ), error(sprintf('(diagn_file) %s does not exist.',diagn_file)), end

% Get some information from the truth_file 
ft = netcdf(truth_file);
t.model      = ft.model(:);
t.num_vars   = ncsize(ft{'StateVariable'}); % determine # of state variables
t.num_copies = ncsize(ft{'copy'}); % determine # of ensemble members
t.num_times  = ncsize(ft{'time'}); % determine # of output times

% Get some information from the diagn_file 
fd = netcdf(diagn_file);
d.model      = fd.model(:);
d.num_vars   = ncsize(fd{'StateVariable'}); % determine # of state variables
d.num_copies = ncsize(fd{'copy'}); % determine # of ensemble members
d.num_times  = ncsize(fd{'time'}); % determine # of output times

% rudimentary bulletproofing
if (strcmp(t.model,d.model) ~= 1)
   error('no No NO ... models must be the same')
end

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(truth_file, 'true state' );
ens_mean_index   = get_copy_index(diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(diagn_file, 'ensemble spread');

% Get the appropriate copies
% (TJH) This function gets them "on the fly" in the plotting.
%truth      = get_state_copy(truth_file, truth_index);
%ens_mean   = get_state_copy(diagn_file, ens_mean_index );
%ens_spread = get_state_copy(diagn_file, ens_spread_index );

% Get some useful plotting arrays
times = getnc(truth_file,'time');

switch lower(t.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3
            ivar = (i - 1)*3 + j;
            disp(sprintf('plotting model %s Variable %d ...',t.model,ivar))
            % Get the truth for this variable
            truth    = get_var_series(truth_file, truth_index, ivar);
            ens_mean = get_var_series(diagn_file, ens_mean_index, ivar);
            subplot(3, 1, j);
            plot(times,truth, 'b',times,ens_mean,'r')
            title(sprintf('model %s Variable %d',t.model,ivar))
            xlabel(sprintf('model time (%d timesteps)',t.num_times))
            s = 'Ensemble Mean';
            legend('True State',s,0)
         end
      end


   case 'lorenz_63'

      disp('lorenz_63 not implemented yet.')

   case 'lorenz_96'
      disp('lorenz_96 not implemented yet.')

   otherwise
end



