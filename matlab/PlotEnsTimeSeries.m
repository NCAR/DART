function PlotEnsTimeSeries(truth_file, diagn_file, state_var_inds )
% PlotEnsTimeSeries: Plots time series of ensemble members, mean and truth.
%
% PlotEnsTimeSeries is intended to be called by 'plot_ens_time_series'
%
% USAGE: PlotEnsTimeSeries(truth_file, diagn_file, state_var_inds);
%
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
% state_var_inds  indices of state variables of interest
%
% Example 1 ( 9var model )
%%--------------------------------------------------------
% truth_file     = 'True_State.nc';
% diagn_file     = 'Posterior_Diag.nc';
% state_var_inds = [ 4 5 6 ];
% PlotEnsTimeSeries(truth_file, diagn_file, state_var_inds);

% TJH Wed Jul  2 09:56:40 MDT 2003

CheckModelCompatibility(truth_file, diagn_file)

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
            truth       = get_var_series(truth_file, truth_index, ivar);
            ens_mean    = get_var_series(diagn_file, ens_mean_index, ivar);
            ens_members = get_ens_series(diagn_file, ivar);  % all members

            subplot(3, 1, j);
               % This is a bit wasteful, but we plot everything once to define
               % the color order for the legend and then again for visibility
               plot(times,      truth,'b','LineWidth',2.0); hold on;
               plot(times,   ens_mean,'r','LineWidth',2.0);
               plot(times,ens_members,'g');
               legend('True State','Ensemble Mean', ...
                      sprintf('Ensemble Members (%d)',d.num_copies-2),0)
               legend boxoff
               plot(times,   truth,'b','LineWidth',2.0); % plot again - on top
               plot(times,ens_mean,'r','LineWidth',2.0); %      again - on top
               title(sprintf('%s Variable %d Ensemble Members of %s',...
                     t.model, ivar, diagn_file), ...
                     'interpreter','none','fontweight','bold')
               xlabel(sprintf('model time (%d timesteps)',t.num_times))
         end
      end

   case 'lorenz_63'

      % Use one figure with three subplots 
      figure(1); clf; iplot = 0;
      for ivar = state_var_inds,

            truth       = get_var_series(truth_file, truth_index, ivar);
            ens_mean    = get_var_series(diagn_file, ens_mean_index, ivar);
            ens_members = get_ens_series(diagn_file, ivar);  % all members

            iplot = iplot + 1;
            subplot(length(state_var_inds), 1, iplot);
               % This is a bit wasteful, but we plot everything once to define
               % the color order for the legend and then again for visibility
               plot(times,      truth,'b','LineWidth',1.0); hold on;
               plot(times,   ens_mean,'r','LineWidth',1.0);
               plot(times,ens_members,'g');
               legend('True State','Ensemble Mean', ...
                      sprintf('Ensemble Members (%d)',d.num_copies-2),0)
               legend boxoff
               plot(times,   truth,'b','LineWidth',1.0); % plot again - on top
               plot(times,ens_mean,'r','LineWidth',1.0); %      again - on top
               title(sprintf('%s Variable %d Ensemble Members of %s',...
                     t.model, ivar, diagn_file), ...
                     'interpreter','none','fontweight','bold')
               xlabel(sprintf('model time (%d timesteps)',t.num_times))
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

      % Use one figure with subplots 
      figure(1); clf; iplot = 0;
      for ivar = state_var_inds,

            truth       = get_var_series(truth_file, truth_index, ivar);
            ens_mean    = get_var_series(diagn_file, ens_mean_index, ivar);
            ens_members = get_ens_series(diagn_file, ivar);  % all members

            iplot = iplot + 1;
            subplot(length(state_var_inds), 1, iplot);
               % This is a bit wasteful, but we plot everything once to define
               % the color order for the legend and then again for visibility
               plot(times,      truth,'b','LineWidth',2); hold on;
               plot(times,   ens_mean,'r','LineWidth',2);
               plot(times,ens_members,'g');
               legend('True State','Ensemble Mean', ...
                      sprintf('Ensemble Members (%d)',d.num_copies-2),0)
               legend boxoff
               plot(times,   truth,'b','LineWidth',2); % plot again - on top
               plot(times,ens_mean,'r','LineWidth',2); %      again - on top
               title(sprintf('%s Variable %d Ensemble Members of %s',...
                     t.model, ivar, diagn_file), ...
                     'interpreter','none','fontweight','bold')
               xlabel(sprintf('model time (%d timesteps)',t.num_times))
      end

   otherwise
      error(sprintf('model %s unknown.',t.model))
end
