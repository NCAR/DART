function PlotBins(truth_file, diagn_file, state_var_inds)

% Plots ensemble rank histograms for several low-order models
%
% Example 1
% truth_file = 'True_State.nc';
% diagn_file = 'Prior_Diag.nc';
% PlotBins(truth_file,diagn_file);

CheckModelCompatibility(truth_file,diagn_file)

% Get the state for the truth
   truth_index = get_copy_index(truth_file,'true state');
ens_mean_index = get_copy_index(diagn_file,'ensemble mean');

true_model = GetAtt(truth_file,'model');

switch lower(true_model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3
            ivar = (i - 1)*3 + j;
            ens   = get_ens_series(diagn_file, ivar );
            truth = get_var_series(truth_file, truth_index, ivar);
            bins  = rank_hist(ens, truth);
            subplot(3, 1, j);
            bar(bins);
            title(sprintf('%s Variable %d for %s', ...
                  true_model,ivar,diagn_file), ...
                  'interpreter','none','fontweight','bold')
         end
      end

   case 'lorenz_63'

      clf; iplot = 0;
      for ivar = state_var_inds,
         iplot = iplot + 1;
         truth = get_var_series(truth_file, truth_index, ivar);
         ens   = get_ens_series(diagn_file, ivar );
         bins  = rank_hist(ens, truth);
         subplot(length(state_var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %d for %s', ...
               true_model,ivar,diagn_file), ...
               'interpreter','none','fontweight','bold')
      end

   case 'lorenz_96'
      % disp('lorenz_96 not implemented yet.')

      clf; iplot = 0;
      for ivar = state_var_inds,
         iplot = iplot + 1;
         truth = get_var_series(truth_file, truth_index, ivar);
         ens   = get_ens_series(diagn_file, ivar );
         bins  = rank_hist(ens, truth);
         subplot(length(state_var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %d for %s', ...
               true_model,ivar,diagn_file), ...
               'interpreter','none','fontweight','bold')
      end


   otherwise
end


%======================================================================
% Subfunctions
%======================================================================

function modelstring = GetAtt(fname,attname)
% Get a global attribute from a netCDF file.

   f = netcdf(fname,'nowrite');   % open with low-level netcdf operators.
   modelstring = f.model(:);   % grab a global attribute
   close(f)
   if isempty(modelstring)
      error(sprintf('NO model in netCDF file %s',fname))
   else
      disp(sprintf('Selected model is %s',modelstring))
   end

