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
   disp(sprintf('%s has model %s ',truth_file,t.model))
   disp(sprintf('%s has model %s ',diagn_file,d.model))
   error('no No NO ... models must be the same')
end
if (t.num_vars ~= d.num_vars)
   disp(sprintf('%s has %d state variables',truth_file,t.num_vars))
   disp(sprintf('%s has %d state variables',diagn_file,d.num_vars))
   error('no No NO ... both files must have same number of state variables.')
end
if (t.num_times ~= d.num_times)
   disp(sprintf('%s has %d timesteps',truth_file,t.num_times))
   disp(sprintf('%s has %d timesteps',diagn_file,d.num_times))
   error('ugh ... both files must have same number of timesteps.')
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
      for i = 1:3,
         figure(i); clf
         for j = 1:3,
            ivar = (i - 1)*3 + j;
            disp(sprintf('plotting model %s Variable %d ...',t.model,ivar))
            % Get the truth for this variable
            truth    = get_var_series(truth_file, truth_index, ivar);
            ens_mean = get_var_series(diagn_file, ens_mean_index, ivar);
            subplot(3, 1, j);
            plot(times,truth, 'b',times,ens_mean,'r')
            title(sprintf('%s Variable %d of %s',t.model,ivar,diagn_file), ...
                  'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',t.num_times))
            s = 'Ensemble Mean';
            legend('True State',s,0)
            legend boxoff
         end
      end

   case 'lorenz_63'

      % Use one figure with three subplots
      clf
      for j = 1:t.num_vars,
            disp(sprintf('plotting model %s Variable %d ...',t.model,j))
            % Get the truth for this variable
            truth    = get_var_series(truth_file, truth_index, j);
            ens_mean = get_var_series(diagn_file, ens_mean_index, j);
            subplot(t.num_vars, 1, j);
            plot(times,truth, 'b',times,ens_mean,'r')
            title(sprintf('%s Variable %d of %s',t.model,j,diagn_file), ...
                  'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',t.num_times))
            s = 'Ensemble Mean';
            legend('True State',s,0)
            legend boxoff
      end
      % as a bonus, plot the mean attractors.
      figure(2); clf
      ts   = get_state_copy(diagn_file,truth_index);
      ens  = get_state_copy(diagn_file,ens_mean_index);
      plot3(  ts(:,1),  ts(:,2),  ts(:,3), 'b', ...
             ens(:,1), ens(:,2), ens(:,3), 'r')
      title(sprintf('%s Attractors for %s and %s', ...
                 t.model, truth_file, diagn_file), ...
                 'interpreter','none','fontweight','bold')
      legend('True State','Ensemble Mean',0)
      legend boxoff
      xlabel('state variable 1')
      ylabel('state variable 2')
      zlabel('state variable 3')

   case 'lorenz_96'
      disp('lorenz_96 not implemented yet.')

   otherwise
end
