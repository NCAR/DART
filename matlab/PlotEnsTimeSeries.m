function PlotEnsTimeSeries(truth_file,diagn_file)
% Plots time series of ensemble members, mean and truth for 9 variable
%
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Posterior_Diag.nc';
% PlotEnsTimeSeries(truth_file,diagn_file);

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

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(truth_file, 'true state' );
ens_mean_index   = get_copy_index(diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(diagn_file, 'ensemble spread');

% Get some useful plotting arrays
times = getnc(truth_file,'time');

switch lower(t.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3

            ivar = (i - 1)*3 + j;
            % disp(sprintf('model %s Variable %d',t.model,ivar))
            truth       = get_var_series(truth_file, truth_index, ivar);
            ens_mean    = get_var_series(diagn_file, ens_mean_index, ivar);
            ens_members = get_ens_series(diagn_file, ivar);  % all members

            subplot(3, 1, j);
               % This is a bit wasteful, but we plot everything once to define
               % the color order for the legend and then again for visibility
               plot(times,      truth,'b','LineWidth',2.0); hold on;
               plot(times,   ens_mean,'r','LineWidth',2.0);
               plot(times,ens_members,'g');
               legend('True State','Ensemble Mean','Ensemble Members',0)
               plot(times,   truth,'b','LineWidth',2.0); % plot again - on top
               plot(times,ens_mean,'r','LineWidth',2.0); %      again - on top
               title(sprintf('model %s Variable %d Ensemble Members (%d)',...
                     t.model, ivar, d.num_copies-2))
               xlabel(sprintf('model time (%d timesteps)',t.num_times))
         end
      end

   case 'lorenz_63'

      disp('lorenz_63 not implemented yet.')

   case 'lorenz_96'
      disp('lorenz_96 not implemented yet.')

   otherwise
end
