function PlotEnsErrSpread(truth_file,diagn_file)
% Plots summary plots of error and spread 
%
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Prior_Diag.nc';
% PlotEnsErrSpread(truth_file,diagn_file)

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
truth      = get_state_copy(truth_file, truth_index);
ens_mean   = get_state_copy(diagn_file, ens_mean_index );
ens_spread = get_state_copy(diagn_file, ens_spread_index );

% Get some useful plotting arrays
times = getnc(truth_file,'time');

switch lower(t.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3

            ivar = (i - 1)*3 + j;

            err         = total_err(ens_mean(:,ivar) , truth(:,ivar));  
            errTotal    = sum(err);
            spreadTotal = sum(ens_spread(:,ivar));
            string1 = ['Ensemble Total Error \Sigma = ' num2str(errTotal)];
            string2 = ['Ensemble Spread \Sigma = ' num2str(spreadTotal)];

            disp(sprintf('%s model Variable %d',t.model,ivar))

            subplot(3, 1, j);
               plot(times,err, 'b', ...
                    times,ens_spread(:, ivar), 'r');
               s1 = sprintf('%s model Var %d Ensemble Error Spread for %s', ...
                            t.model,ivar,diagn_file);
               title(s1,'interpreter','none','fontweight','bold')
               legend(string1,string2,0)
               legend boxoff
               xlabel(sprintf('model time (%d timesteps)',t.num_times))
         end
      end

   case 'lorenz_63'

      % disp('lorenz_63 not implemented yet.')
      % Use one figure with three subplots
      clf
      for j = 1:t.num_vars

            err         = total_err(ens_mean(:,j) , truth(:,j));  
            errTotal    = sum(err);
            spreadTotal = sum(ens_spread(:,j));
            string1 = ['Ensemble Total Error \Sigma = ' num2str(errTotal)];
            string2 = ['Ensemble Spread \Sigma = ' num2str(spreadTotal)];

            disp(sprintf('model %s Variable %d',t.model,j))

            subplot(t.num_vars, 1, j);
               plot(times,err, 'b', ...
                    times,ens_spread(:, j), 'r');
               s1 = sprintf('%s model Var %d Ensemble Error Spread for %s', ...
                            t.model,j,diagn_file);
               title(s1,'interpreter','none','fontweight','bold')
               legend(string1,string2,0)
               legend boxoff
               xlabel(sprintf('model time (%d timesteps)',t.num_times))
      end

   case 'lorenz_96'
      disp('lorenz_96 not implemented yet.')

   otherwise
end
