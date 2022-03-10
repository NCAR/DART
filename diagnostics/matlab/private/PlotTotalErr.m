function PlotTotalErr( pinfo )
%% PlotTotalErr Plots summary plots of global error and spread
%
% The error is a measure of the distance between the
% 'ensemble mean' and the 'true_state'. The distance of the
% spread is also plotted.
%
% PlotTotalErr is intended to be called by 'plot_total_err'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotTotalErr( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
%
% Example 1   (Lorenz_63, Lorenz_96, 9var ...)
%%--------------------------------------------------------
% pinfo.truth_file = 'true_state.nc';
% pinfo.diagn_file = 'postassim.nc';
% PlotTotalErr( pinfo )

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

figdata = setfigure('landscape');

switch lower(pinfo.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04','ikeda', 'null'}

      %% Get the appropriate netcdf variables

      truth  = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', 'state', ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'tindex1',pinfo.truth_time(1), ...
                   'tcount',pinfo.truth_time(2));
      ens    = get_hyperslab('fname',pinfo.diagn_file, ...
                   'varname','state_mean', ...
                   'permute', 'T', ...
                   'tindex1',pinfo.diagn_time(1), ...
                   'tcount',pinfo.diagn_time(2));
      spread = get_hyperslab('fname',pinfo.diagn_file, ...
                   'varname','state_sd', ...
                   'permute', 'T', ...
                   'tindex1',pinfo.diagn_time(1), ...
                   'tcount',pinfo.diagn_time(2));

      [~,num_vars] = size(spread);

      % Also need to compute the spread; zero truth for this and
      % compute distance from 0
      err        = squeeze(total_err(truth, ens));
      err_spread = squeeze(total_err(zeros(size(spread)), spread));
      errTotal   = sum(err)/pinfo.time_series_length;
      spreadTotal= sum(err_spread)/pinfo.time_series_length;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

      clf;
      plot(pinfo.time,err, 'b', pinfo.time,err_spread, 'r','LineWidth',figdata.linewidth);
      set(gca,'FontSize',figdata.fontsize)
      legend(string1,string2,'Location','NorthEast','FontSize',figdata.fontsize)
      legend boxoff
      s1 = sprintf('%s Total Error over all %d variables', pinfo.model, num_vars);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
      ylabel('Total Error')

   case 'lorenz_96_2scale'

      %% Simply going to append X,Y together and treat as above.
      % Get the appropriate netcdf variables

      tim    = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', 'X', ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'tindex1', pinfo.truth_time(1), ...
                   'tcount', pinfo.truth_time(2));
      tom    = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', 'Y', ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'tindex1', pinfo.truth_time(1), ...
                   'tcount', pinfo.truth_time(2));

      truth  = [tim tom];

      tim    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', 'X_mean', ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));
      tom    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', 'Y_mean', ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

      ens    = [tim tom];

      tim    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', 'X_sd', ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));
      tom    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', 'Y_sd', ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

      spread = [tim tom]; clear tim tom
      [~,num_vars] = size(spread);

      % Also need to compute the spread; zero truth for this and
      % compute distance from 0
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/pinfo.time_series_length;
      spreadTotal= sum(err_spread)/pinfo.time_series_length;
      string1 = sprintf('time-mean Ensemble Mean Total Error = %f',errTotal);
      string2 = sprintf('time-mean Ensemble Spread = %f',spreadTotal);

      clf;
      plot(pinfo.time,err, 'b', pinfo.time ,err_spread, 'r','LineWidth',figdata.linewidth);
      set(gca,'FontSize',figdata.fontsize)
      legend(string1,string2,'Location','NorthEast','FontSize',figdata.fontsize)
      legend boxoff
      s1 = sprintf('%s Total Error over all %d variables', pinfo.model, num_vars);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
      ylabel('Total Error')

   case 'forced_lorenz_96'

      %% This model has the state variables replicated, so there is a difference
      % between num_state_vars and the length of the state variable.

      num_model_vars = ncreadatt(pinfo.truth_file, '/', 'model_num_state_vars');

      % Get the appropriate netcdf variables

      Whole_truth  = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', 'state', ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'tindex1', pinfo.truth_time(1), ...
                   'tcount', pinfo.truth_time(2));

      Whole_ens    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', 'state_mean', ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

      Whole_spread = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', 'state_sd', ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

      [~,num_vars] = size(Whole_spread);

      %--------------------------------------------------------------------------
      % Treat the traditional state variable independent of the forcing variables
      %--------------------------------------------------------------------------

      ind1 = 1;                 % ASSUMPTION: traditional state is first
      indN = num_model_vars;

      truth  = Whole_truth(  :, ind1:indN );
      ens    = Whole_ens(    :, ind1:indN );
      spread = Whole_spread( :, ind1:indN );

      % Compute the spread; zero truth for this and compute distance from 0

      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/pinfo.time_series_length;
      spreadTotal= sum(err_spread)/pinfo.time_series_length;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

      clf; subplot(2,1,1);

      plot(pinfo.time,err, 'b', pinfo.time,err_spread, 'r','LineWidth',figdata.linewidth);
      set(gca,'FontSize',figdata.fontsize)
      legend(string1,string2,'Location','NorthEast','FontSize',figdata.fontsize)
      legend boxoff
      s1 = sprintf('%s Total Error over statevars %d to %d', pinfo.model, ind1, indN);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
      ylabel('Total Error')

      %--------------------------------------------------------------------------
      % Now for the forcing
      %--------------------------------------------------------------------------

      ind1 = num_model_vars + 1;
      indN = num_vars;

      truth  = Whole_truth(  :, ind1:indN );
      ens    = Whole_ens(    :, ind1:indN );
      spread = Whole_spread( :, ind1:indN );

      % Compute the spread; zero truth for this and compute distance from 0

      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/pinfo.time_series_length;
      spreadTotal= sum(err_spread)/pinfo.time_series_length;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

      subplot(2,1,2)

      plot(pinfo.time,err, 'b', pinfo.time,err_spread, 'r','LineWidth',figdata.linewidth);
      set(gca,'FontSize',figdata.fontsize)
      legend(string1,string2,'Location','NorthEast','FontSize',figdata.fontsize)
      legend boxoff
      s1 = sprintf('%s Total Error over statevars %d to %d', pinfo.model, ind1, indN);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
      ylabel('Total Error')

   case {'simple_advection'}

      varlist = {'concentration','source','wind','mean_source','source_phase'};

      for ivar = 1:length(varlist)

         % Get the appropriate netcdf variables

         truth  = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', varlist{ivar}, ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'tindex1', pinfo.truth_time(1), ...
                   'tcount', pinfo.truth_time(2));

         varname = sprintf('%s_mean',varlist{ivar});
         ens    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', varname, ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

         varname = sprintf('%s_sd',varlist{ivar});
         spread = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', varname, ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));
               
         [~,num_vars] = size(spread);

         % Also need to compute the spread; zero truth for this and
         % compute distance from 0
         err        = total_err(truth, ens);
         err_spread = total_err(zeros(size(spread)), spread);
         errTotal   = sum(err)/pinfo.time_series_length;
         spreadTotal= sum(err_spread)/pinfo.time_series_length;
         string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
         string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

         figure(ivar); clf(ivar);
         plot(pinfo.time,err, 'b', pinfo.time,err_spread, 'r','LineWidth',figdata.linewidth);
         set(gca,'FontSize',figdata.fontsize)
         legend(string1,string2,'Location','NorthEast','FontSize',figdata.fontsize)
         legend boxoff
         string1 = sprintf('%s Total Error over all %d variables', pinfo.model, num_vars);
         string2 = sprintf('''%s'' %s', varlist{ivar}, pinfo.diagn_file);
         title({string1,string2},'interpreter','none','fontweight','bold')
         xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
         ylabel('Total Error')
      end
      
      
   case {'lorenz_96_tracer_advection'}

      varlist = {'state_variable','tracer_concentration','source'};

      for ivar = 1:length(varlist)

         % Get the appropriate netcdf variables

         truth  = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', varlist{ivar}, ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'tindex1', pinfo.truth_time(1), ...
                   'tcount', pinfo.truth_time(2));

         varname = sprintf('%s_mean',varlist{ivar});
         ens    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', varname, ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

         varname = sprintf('%s_sd',varlist{ivar});
         spread = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', varname, ...
                   'permute', 'T', ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));
               
         [~,num_vars] = size(spread);

         % Also need to compute the spread; zero truth for this and
         % compute distance from 0
         err        = total_err(truth, ens);
         err_spread = total_err(zeros(size(spread)), spread);
         errTotal   = sum(err)/pinfo.time_series_length;
         spreadTotal= sum(err_spread)/pinfo.time_series_length;
         string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
         string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

         figure(ivar); clf(ivar);
         plot(pinfo.time,err, 'b', pinfo.time,err_spread, 'r','LineWidth',figdata.linewidth);
         set(gca,'FontSize',figdata.fontsize)
         legend(string1,string2,'Location','NorthEast','FontSize',figdata.fontsize)
         legend boxoff
         string1 = sprintf('%s Total Error over all %d variables', pinfo.model, num_vars);
         string2 = sprintf('''%s'' %s', varlist{ivar}, pinfo.diagn_file);
         title({string1,string2},'interpreter','none','fontweight','bold')
         xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
         ylabel('Total Error')
      end

   case 'fms_bgrid'
      %% GFDL bgrid model

      BgridTotalError( pinfo )

   case 'pe2lyr'
      error('not supported yet')
      %% primitive equation 2 layer model

      Pe2lyrTotalError( pinfo )

   case 'pbl_1d'
      error('not supported yet')
      %% planetary boundary layer column model

      PBL1DTotalError( pinfo )

   case 'mitgcm_ocean'
      error('not supported yet')
      %% MIT general circulation ocean model

      MITGCMOceanTotalError( pinfo )

   case 'cam'
      error('not supported yet')
      %% community (global) general circulation atmosphere model

      CAMTotalError( pinfo )

   case 'wrf'
      error('not supported yet')
      %% weather and research forecasting model

      WRFTotalError( pinfo )

   case 'mpas_atm'
      error('not supported yet')
      %% unstructured grid atmosphere model

      MPAS_ATMTotalError( pinfo )

   case 'pop'
      error('not supported yet')
      %% parallel ocean program

      POPTotalError( pinfo )

   case 'sqg'
      error('not supported yet')

      SqgTotalError( pinfo )

   otherwise

      fprintf('unsupported model %s -- doing nothing\n',pinfo.model)

end

%=======================================================================
% helper functions
%=======================================================================

function PBL1DTotalError ( pinfo )

%% Get some standard plotting arrays

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

error('this block is deprecated, use get_hyperslab() instead.')

count = [pinfo.time_series_length 1];

start  = [pinfo.truth_time(1)-1  truth_index-1];
truth  = nc_varget(pinfo.truth_file,'U', start, count);

start  = [pinfo.diagn_time(1)-1  ens_mean_index-1];
ens    = nc_varget(pinfo.diagn_file,'U', start, count);

start  = [pinfo.diagn_time(1)-1  ens_spread_index-1];
spread = nc_varget(pinfo.diagn_file,'U', start, count);

err        = total_err(              truth,    ens);
err_spread = total_err(zeros(size(spread)), spread);

y_error  = squeeze(mean(err,2));           % mean over all levels
y_spread = squeeze(mean(err_spread,2));    % mean over all levels

plot(pinfo.time,y_error,'r-',pinfo.time,y_spread,'g-')
title('PBL_1d mean error of U over time ... all levels.')
xdates(pinfo.time)
ylabel('mean (all levels) total error')
axis([-Inf Inf 0 Inf])




function MITGCMOceanTotalError( pinfo )
%% -------------------------------------------------------------------
% netcdf has only prognostic variables.
% We are going to plot the total error (over a horizontal slice)
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

tstart = pinfo.truth_time(1);
tcount = pinfo.truth_time(2);
dstart = pinfo.diagn_time(1);
dcount = pinfo.diagn_time(2);

error('this block is deprecated, use get_hyperslab() instead.')
% Since the models are "compatible", get the info from either one.
XC        = nc_varget(pinfo.truth_file, 'XC'); num_XC  = length(XC );
YC        = nc_varget(pinfo.truth_file, 'YC'); num_YC  = length(YC );
XG        = nc_varget(pinfo.truth_file, 'XG');
YG        = nc_varget(pinfo.truth_file, 'YG');
ZC        = nc_varget(pinfo.truth_file, 'ZC'); num_ZC  = length(ZC);
num_times = pinfo.time_series_length;

% Initialize storage for error averaging
rms          = zeros(num_times, pinfo.num_state_vars, num_ZC);
spread_final = zeros(num_times, pinfo.num_state_vars, num_ZC);

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = reshape(SphereWeights(XC, YC),1,num_YC*num_XC);   % Grid Centers
vwts = reshape(SphereWeights(XG, YG),1,num_YC*num_XC);   % Grid edGes

%----------------------------------------------------------------------
% temperature ...  num_times x num_levels x num_lats x num_lons
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_ZC,     % Loop through all levels
for ivar=1:pinfo.num_state_vars-1,  % knowing ssh is last

   fprintf('Processing level %d of %d ...\n',ilevel,num_ZC)
   %-------------------------------------------------------------------
   % all vars organized    num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------

   truth  = GetLevel(pinfo.truth_file, pinfo.vars{ivar},      truth_index, ilevel, tstart, tcount);
   ens    = GetLevel(pinfo.diagn_file, pinfo.vars{ivar},   ens_mean_index, ilevel, dstart, dcount);
   spread = GetLevel(pinfo.diagn_file, pinfo.vars{ivar}, ens_spread_index, ilevel, dstart, dcount);

   landmask = find(isfinite(truth(1,:))); % presume all members have same mask
   vweights = vwts(landmask);
   tweights = twts(landmask);

   mytruth  =  truth(:,landmask);
   myens    =    ens(:,landmask);
   myspread = spread(:,landmask);

   switch lower(pinfo.vars{ivar})
      case {'v'}
         err        = total_err(              mytruth,    myens, vweights);
         err_spread = total_err(zeros(size(myspread)), myspread, vweights);
      otherwise
         err        = total_err(              mytruth,    myens, tweights);
         err_spread = total_err(zeros(size(myspread)), myspread, tweights);
   end

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;

end
end % End of level loop

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% surface pressure has only one level.
% Get2D   returns a   num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

disp('Processing sea surface height ...')

ivar   = find(strcmp(pinfo.vars,'SSH'));
ilevel = 1;

truth      = Get2D(pinfo.truth_file, 'SSH',      truth_index, tstart, tcount);
ens        = Get2D(pinfo.diagn_file, 'SSH',   ens_mean_index, dstart, dcount);
spread     = Get2D(pinfo.diagn_file, 'SSH', ens_spread_index, dstart, dcount);

landmask   = find(isfinite(truth(1,:))); % presume all members have same mask
tweights   = twts(landmask);

mytruth    =  truth(:,landmask);
myens      =    ens(:,landmask);
myspread   = spread(:,landmask);

err        = total_err(              mytruth,    myens, tweights );
err_spread = total_err(zeros(size(myspread)), myspread, tweights );

         rms(:,ivar,ilevel) = err;         % spatial mean
spread_final(:,ivar,ilevel) = err_spread;  % spatial mean

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% Each variable in its own figure window
%----------------------------------------------------------------------
for ivar=1:pinfo.num_state_vars,

   figure(ivar); clf;
      varunits = ncreadatt(pinfo.truth_file,pinfo.vars{ivar},'units');

      s1 = sprintf('%s %s Ensemble Mean', pinfo.model,pinfo.vars{ivar});

      switch lower(pinfo.vars{ivar})
         case {'ssh'}
            plot(pinfo.time,          rms(:, ivar, 1), '-'); hold on;
            plot(pinfo.time, spread_final(:, ivar, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean Total Error = %f', mean(         rms(:, ivar, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
         otherwise
            plot(pinfo.time, squeeze(         rms(:, ivar, :)),'-'); hold on;
            plot(pinfo.time, squeeze(spread_final(:, ivar, :)),'--');

            for i = 1:num_ZC,
               s{i       } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
               s{i+num_ZC} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
            end
      end

      set(gca,'FontSize',figdata.fontsize)
      h = legend(s); legend(h,'boxoff','FontSize',figdata.fontsize)
      grid on;
      xdates(pinfo.time)
      ylabel(sprintf('global-area-weighted rmse (%s)',varunits))
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end



function Pe2lyrTotalError( pinfo )
%% -------------------------------------------------------------------
% netcdf has no state vector, it has prognostic variables.
% Plot the area-weighted total error for each variable and annotate a
% time-mean value.
%---------------------------------------------------------------------

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
rads    = zeros(numel(pinfo.lats),1);
rads(:) = pi*pinfo.lats/180.0;
wts     = cos(rads) * ones(1,numel(pinfo.lons));
wts     = wts ./ sum(wts(:));
weights = wts(:);

%----------------------------------------------------------------------
%
%----------------------------------------------------------------------

for ivar=1:pinfo.num_state_vars,

   fprintf('Processing %s ...\n',pinfo.vars{ivar})
   varunits = ncreadatt(pinfo.truth_file, pinfo.vars{ivar}, 'units');
   rmse     = zeros(pinfo.time_series_length,1);
   sprd     = zeros(pinfo.time_series_length,1);

   %%------------------------------------------------------------------
   % all vars organized     num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   for itime = 1:pinfo.time_series_length,

      truth  = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'copyindex', truth_index, ...
                   'timeindex', pinfo.truth_time(1)+itime-1);

      ens    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'copyindex', ens_mean_index, ...
                   'timeindex', pinfo.diagn_time(1)+itime-1);

      spread = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'copyindex', ens_spread_index, ...
                   'timeindex', pinfo.diagn_time(1)+itime-1);

      sqerr       = (truth - ens).^2;
      sqsprd      =     spread   .^2;
      flatsqerr   = squeeze(mean(sqerr,1));         % average over levels
      flatsqsprd  = squeeze(mean(sqsprd,1));        % average over levels
      ms_err      = sum(flatsqerr(:)  .* weights);
      ms_spread   = sum(flatsqsprd(:) .* weights);
      rmse(itime) = sqrt(ms_err);
      sprd(itime) = sqrt(ms_spread);
   end

   clear truth ens spread sqerr sqsprd flatsqerr flatsqsprd ms_err ms_spread

   %-------------------------------------------------------------------
   % Each variable in its own figure window
   %-------------------------------------------------------------------

   figure(ivar); clf;

   plot(pinfo.time, rmse,'-', pinfo.time, sprd,'--');

   s1   = sprintf('%s Ensemble Mean for %s', pinfo.model,pinfo.vars{ivar});
   s{1} = sprintf('time-mean Ensemble Mean Total Error  = %f', mean(rmse));
   s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(sprd));

   set(gca,'FontSize',figdata.fontsize)
   h = legend(s); legend(h,'boxoff','FontSize',figdata.fontsize)
   grid on;
   xdates(pinfo.time)
   ylabel(sprintf('global-area-weighted rmse (%s)',varunits))
   title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end


function MPAS_ATMTotalError( pinfo )
%% -------------------------------------------------------------------
% Assume netcdf has only prognostic variables.
% We are going to plot a time series of the total area-weighted error
% for each variable (all levels combined and considered equal).
%---------------------------------------------------------------------

num_times  = pinfo.time_series_length;

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Generate the area-based weights. Since we've got a hybrid vertical
% coordinate system, this is a simplification that is not perfect.
weights = pinfo.area ./ sum(pinfo.area(:));

%----------------------------------------------------------------------
% get_hyperslab returns a  num_times x num_gridcells ... usually.
% If time is a singleton dimension, it returns    num_gridcells x 1
%----------------------------------------------------------------------

for ivar=1:pinfo.num_state_vars,

   varunits = ncreadatt(pinfo.truth_file, pinfo.vars{ivar}, 'units');
   s1 = sprintf('%s %s Ensemble Mean', pinfo.model, pinfo.vars{ivar});
   fprintf('%s ...\n',s1)

   % Initialize storage for error averaging
   num_levels   = FindNLevels(   pinfo.truth_file, pinfo.vars{ivar});
   num_locations= FindNLocations(pinfo.truth_file, pinfo.vars{ivar});
   sqerr        = zeros(num_times, num_locations);
   sqspread     = zeros(num_times, num_locations);

   if (length(weights) ~= num_locations)
      error('%s is not a cell-based quantity. Cannot apply weights.\n',pinfo.vars{ivar})
   end

   %% Loop over all levels for this variable
   for ilevel = 1:num_levels,

      %-------------------------------------------------------------------
      % want all vars organized num_times x num_gridcells
      %-------------------------------------------------------------------

      truth  = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'levelindex', ilevel, ...
                   'copyindex', truth_index, ...
                   'tindex1', pinfo.truth_time(1), ...
                   'tcount', pinfo.truth_time(2));

      ens    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'levelindex', ilevel, ...
                   'copyindex', ens_mean_index, ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

      spread = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'levelindex', ilevel, ...
                   'copyindex', ens_spread_index, ...
                   'tindex1',  pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

      % if there is one timestep, the vars get returned transposed.
      if (pinfo.truth_time(2) == 1)
          truth  = truth';
          ens    = ens';
          spread = spread';
      end

      sqerr    = sqerr    + (truth - ens).^2;
      sqspread = sqspread + spread .^2;

   end

   clear truth ens spread

   sqerr     = sqerr    / num_levels; % average of all levels
   sqspread  = sqspread / num_levels;
   rmsetot   = sqrt( sqerr    * weights );
   spreadtot = sqrt( sqspread * weights );

   %% Plot the time-evolution for all levels on one figure.

   figure(ivar); clf;

   plot(pinfo.time, rmsetot,'-',pinfo.time, spreadtot,'--');

   s = {sprintf('time-mean Ensemble Mean Total Error = %f', mean(rmsetot)), ...
        sprintf('time-mean Ensemble Spread = %f', mean(spreadtot))};

   set(gca,'FontSize',figdata.fontsize)
   h = legend(s); legend(h,'boxoff','FontSize',figdata.fontsize)
   grid on;
   xdates(pinfo.time)
   ylabel(sprintf('global-area-weighted rmse (%s)',varunits))
   title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end


function SqgTotalError( pinfo )
%% -------------------------------------------------------------------
% netcdf has no state vector, it has prognostic variables.
% We are going to plot the total error (over a horizontal slice)
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

model     = ncreadatt(pinfo.truth_file, '/', 'model');
timeunits = ncreadatt(pinfo.truth_file, 'time', 'units');

tstart = pinfo.truth_time(1);
tcount = pinfo.truth_time(2) - pinfo.truth_time(1) + 1;
dstart = pinfo.diagn_time(1);
dcount = pinfo.diagn_time(2) - pinfo.diagn_time(1) + 1;

error('this block is deprecated, use get_hyperslab() instead.')
% Since the models are "compatible", get the info from either one.
tlons    = nc_varget(pinfo.truth_file, 'lon' ); num_tlons  = length(tlons );
tlats    = nc_varget(pinfo.truth_file, 'lat' ); num_tlats  = length(tlats );
levels   = nc_varget(pinfo.truth_file, 'lev' ); num_levels = length(levels);
ens_mems = nc_varget(pinfo.diagn_file, 'copy'); ens_size   = length(ens_mems);

times         = nc_varget(pinfo.truth_file, 'time', tstart-1, tcount);
num_times     = length(times );

% Try to coordinate "time" ... a poor attempt, needs refining
ens_times     = nc_varget(pinfo.diagn_file, 'time', dstart-1, dcount);
num_ens_times = length(ens_times);
if num_ens_times < num_times
   times     =     ens_times;
   num_times = num_ens_times;
end

% Initialize storage for error averaging
rms          = zeros(num_times, pinfo.num_state_vars, num_levels);
spread_final = zeros(num_times, pinfo.num_state_vars, num_levels);

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = reshape(SphereWeights(tlons, tlats),1,num_tlats*num_tlons);

%----------------------------------------------------------------------
% Trying not to assume we can get the whole 3D array at once.
% total_err needs a   num_times-x-num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_levels,
for ivar=1:pinfo.num_state_vars,

   fprintf('Processing level %d of %d ...\n',ilevel,num_levels)

   truth3D  = get_hyperslab('fname', pinfo.truth_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'squeeze', 'T', ...
                   'levelindex', ilevel, ...
                   'copyindex', truth_index, ...
                   'tindex1', pinfo.truth_time(1), ...
                   'tcount', pinfo.truth_time(2));

   ens3D    = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'levelindex', ilevel, ...
                   'copyindex', ens_mean_index, ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

   spread3D = get_hyperslab('fname', pinfo.diagn_file, ...
                   'varname', pinfo.vars{ivar}, ...
                   'permute', 'T', ...
                   'levelindex', ilevel, ...
                   'copyindex', ens_spread_index, ...
                   'tindex1', pinfo.diagn_time(1), ...
                   'tcount', pinfo.diagn_time(2));

   truth  = reshape( truth3D, num_times, num_tlats*num_tlons); clear truth3D
   ens    = reshape(   ens3D, num_times, num_tlats*num_tlons); clear ens3D
   spread = reshape(spread3D, num_times, num_tlats*num_tlons); clear spread3D

   err        = total_err(              truth,    ens, twts);
   err_spread = total_err(zeros(size(spread)), spread, twts);

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;

end
end

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% Each variable in its own figure window
%----------------------------------------------------------------------
for ivar=1:pinfo.num_state_vars,

   figure(ivar); clf;

      varunits = ncreadatt(pinfo.truth_file, pinfo.vars{ivar}, 'units');

      s1 = sprintf('%s %s Ensemble Mean for %s', model,pinfo.vars{ivar});

      switch lower(pinfo.vars{ivar})
         case {'ps'}
            h1 = plot(times,          rms(:, ivar, 1), '-'); hold on;
            h2 = plot(times, spread_final(:, ivar, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, ivar, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
         otherwise
            h1 = plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
            h2 = plot(times, squeeze(spread_final(:, ivar, :)),'--');

            for i = 1:num_levels,
               s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
               s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
            end
      end

      set(gca,'FontSize',figdata.fontsize)
      %h = legend([h1 h2],s); legend(h,'boxoff')
      h = legend(s); legend(h,'boxoff','FontSize',figdata.fontsize)
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end


function POPTotalError( pinfo )
%% -------------------------------------------------------------------
% Lets face it - the total error is completely ill-defined for large
% models. Just issue a statement and abort.
%---------------------------------------------------------------------
disp(' ')
disp('There is no meaningful way to calculate total error for POP.')
disp('Should each variable be considered separately?')
disp('Should we normalize each variable - What about each level?')
disp('The levels do not have the same number of ocean gridcells  ...')
disp('Consequently, we are doing nothing. Sorry.')
disp(' ')


function slice = Get2D(fname, varname, copyindex, tstartind, tcount )
%% this gets a bit funky if tcount == 1;
% it automatically squeezes the singleton dimension,
% leaving the dims in the wrong slots.  copyindex will always
% squeeze out, so try to figure out if you get a 3d or 2d return
% back from getnc() and then deal with time explicitly.
% same goes below in GetLevel().

disp('Get2D() is deprecated, use get_hyperslab() instead.')

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tcount;
         break
      otherwise
   end
end
ted = nc_varget(fname, varname, start, count);
n   = ndims(ted);
if (n == 2)
   [ny,nx] = size(ted);
   nt = 1;
else
   [nt,ny,nx] = size(ted);
end
slice      = reshape(ted,[nt ny*nx]);




function slice = GetLevel(fname, varname, copyindex, ilevel, tstartind, tcount)
%%

disp('GetLevel() is deprecated, use get_hyperslab() instead.')

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.levelindex = ilevel;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tcount;
         break
      otherwise
   end
end

ted = nc_varget(fname, varname, start, count);
n   = ndims(ted);
if (n == 2)
   [ny,nx] = size(ted);
   nt = 1;
else
   [nt,ny,nx] = size(ted);
end
slice      = reshape(ted,[nt ny*nx]);



function wts = SphereWeights(lons, lats)
%% SphereWeights creates weights based on area ...
%
% lons,lats must be 1D arrays (in degrees)

nlons = length(lons);
nlats = length(lats);

if ( numel(lons) ~= nlons )
   disp('longitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end
if ( numel(lats) ~= nlats )
   disp('latitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end

rads    = zeros(nlats,1);               % Ensure lats is a column vector,
rads(:) = pi*lats/180.0;                % and convert to radians.
wts     = cos( rads ) * ones(1,nlons);  % Results in a [nlat-x-nlon] matrix.
wts     = wts ./ sum(wts(:));           % Normalize to unity.



function nlevels = FindNLevels(fname, varname)
%%
varinfo = nc_getvarinfo(fname,varname);
nlevels = 1; % implicitly surface-only variable

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'nvertlevels','nvertlevelsp1'}
         nlevels = varinfo.Size(i);
         break
      otherwise
   end
end



function nlocations = FindNLocations(fname, varname)
%%
varinfo    = nc_getvarinfo(fname,varname);
nlocations = [];

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'ncells','nedges'}
         nlocations = varinfo.Size(i);
         break
      otherwise
   end
end



function xdates(dates)
if (length(dates) < 5)
   set(gca,'XTick',dates);
   datetick('x',31,'keepticks','keeplimits');
   xlabel('Model date (YYYY-MM-DD HH:MM:SS)')
else
   datetick('x','mm.dd.HH','keeplimits'); % 'mm/dd'
   monstr = datestr(dates(1),31);
   xlabel(sprintf('month.day.HH - %s start',monstr))
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
