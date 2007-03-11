function PlotSawtooth( pinfo )
% PlotSawtooth: Plots time series of ensemble members, mean and truth.
%
%
% PlotSawtooth is intended to be called by 'plot_sawtooth'
% The only input argument is a structure with model-dependent
% components.
%
% Both the prior and posterior estimates are plotted as a single
% trajectory. If there is no change in the model state, this should
% appear as a series of steps. This necessitates plotting the 'posterior'
% first ... think about it ...
%
% USAGE: PlotSawtooth( pinfo );
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% prior_file      name of forecast netCDF DART file with desired copy'
% posterior_file  name of nudged netCDF DART file with desired copy'
% var             name of netCDF variable of interest
% var_inds        indices of variables of interest
%
% Example 1 ( forced_lorenz_96 model )
%%--------------------------------------------------------
% pinfo.truth_file     = 'True_State.nc';
% pinfo.prior_file     = 'Prior_Diag.nc';
% pinfo.posterior_file = 'Posterior_Diag.nc';
% pinfo.var            = 'state';
% pinfo.var_inds       = [ 23 36 42 ]; 
% PlotSawtooth( pinfo );
%
% Example 2 ( fms_bgrid_model )
%%--------------------------------------------------------
% pinfo.truth_file     = 'True_State.nc';
% pinfo.prior_file     = 'Prior_Diag.nc';
% pinfo.posterior_file = 'Posterior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotSawtooth( pinfo )

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

CheckModelCompatibility(pinfo.prior_file, pinfo.posterior_file)

% Get some information from the truth_file, if it exists.
if ( exist(pinfo.truth_file,'file') == 2 ) 
   CheckModelCompatibility(pinfo.truth_file, pinfo.posterior_file)

   ft = netcdf(pinfo.truth_file);
   truth.model      = ft.model(:);
   truth.num_times  = length(ft('time'));          % determine # of output times
   close(ft);

   truth.truth_index  = get_copy_index(pinfo.truth_file, 'true state' );
else
   truth = [];
end

% Get some information from the prior_file 
fpr = netcdf(pinfo.prior_file);
prior.model      = fpr.model(:);
prior.num_times  = length(fpr('time'));          % determine # of output times
prior.times      = getnc(pinfo.prior_file,'time');
close(fpr);

% Get some information from the posterior_file 
fpo = netcdf(pinfo.posterior_file);
post.model      = fpo.model(:);
post.num_times  = length(fpo('time'));          % determine # of output times
post.timeunits  = fpo{'time'}.units(:);
post.calendar   = fpo{'time'}.calendar(:);
close(fpr);

% Get the indices for the true state and ensemble mean
% The metadata is queried to determine which "copy" is appropriate.
prior.ens_mean_index = get_copy_index(pinfo.prior_file,     'ensemble mean');
 post.ens_mean_index = get_copy_index(pinfo.posterior_file, 'ensemble mean');

% Get some useful plotting arrays
ens_mean   = zeros(2,prior.num_times);   % trajectory of ensemble mean
ens_member = zeros(2,prior.num_times);   % trajectory of ensemble member
x          = zeros(2,prior.num_times);   % xaxis 
x(1,:)     = prior.times;
x(2,:)     = prior.times;
pinfo.xax  = x(:);
metadata   = getnc(pinfo.prior_file,'CopyMetaData'); % get all the metadata

if isfield(pinfo,'var_inds')
   PlotGivenIndices( pinfo, prior, post, truth, metadata);
elseif isfield(pinfo,'var_names')
   PlotGivenVariable(pinfo, prior, post, truth, metadata);
end

%------------------------------------------------------------------------------
% Plot given an index into the state-space vector
%------------------------------------------------------------------------------

function PlotGivenIndices(pinfo, prior, post, truth, metadata)

   nfigs = length(pinfo.var_inds);  % each variable gets its own figure
   iplot = 0;

   for ivar = pinfo.var_inds,

      % Get the data from the netcdf files

      po_ens_mean = get_var_series(pinfo.posterior_file, pinfo.var, ...
                                   post.ens_mean_index, ivar);
      pr_ens_mean = get_var_series(pinfo.prior_file, pinfo.var, ...
                                   prior.ens_mean_index, ivar);

      % Now we paste them together in a clever way to show 
      % the effect of the assimilation
      ens_mean(  1,:) = pr_ens_mean;
      ens_mean(  2,:) = po_ens_mean;
      a               = ens_mean(:);

      % Plot the true trajectory if it exists; the ens mean; annotate

      iplot = iplot + 1;
      figure(iplot); clf; 

      if (isfield(pinfo,'truth_file'))
         true_trajectory = get_var_series(pinfo.truth_file, pinfo.var, truth.truth_index, ivar);
	 h = plot(prior.times, true_trajectory, 'k-','linewidth',1.0); hold on;
	 h = plot(pinfo.xax, a, 'k-','linewidth',2.0);
	 legend('truth','ensemble mean')
      else 
	 h = plot(pinfo.xax, a, 'k-','linewidth',2.0); hold on;
	 legend('ensemble mean')
      end

      ylabel(sprintf('''%s'' index %d', pinfo.var, ivar ))
      xlabel(sprintf('model time (%d timesteps)', prior.num_times ))
      title(sprintf('%s Trajectories',pinfo.model), ...  
                     'interpreter','none','fontweight','bold')

      % Now check to see if we are overlaying any individual ensemble members.
      % if pinfo.copyindices = [], nothing happens. 

      ens_colors = get(gca,'ColorOrder');   % trying to cycle through colors
      ncolors = size(ens_colors,1) - 1;     % last one is black, already used.

      nmem = 0;
      for imem = pinfo.copyindices
         nmem = nmem + 1;

	 %str1 = sprintf('ensemble member %d',imem);
	 str1 = deblank(metadata(imem,:));
         copy_index = get_copy_index(pinfo.prior_file, str1);
         po_series  = get_var_series(pinfo.posterior_file, pinfo.var, copy_index, ivar);
         pr_series  = get_var_series(pinfo.prior_file,     pinfo.var, copy_index, ivar);

	 ens_member(1,:) = pr_series;
	 ens_member(2,:) = po_series;
	 b               = ens_member(:);

   	 hold on; 
	 memcolor = 1 + mod(nmem-1,ncolors); % cycles through colors [1,6]
	 h = plot(pinfo.xax, b,'linewidth',0.5,'Color',ens_colors(memcolor,:));

	 [legh, objh, outh, outm] = legend;
	 nlines = length(outm);
         outm{nlines+1} = str1;
	 [legh, objh, outh, outm] = legend([outh; h],outm,0);
      end
      legend boxoff

   end

%------------------------------------------------------------------------------
% Plot given an index into the state-space vector
%------------------------------------------------------------------------------

function PlotGivenVariable(pinfo, prior, post, truth, metadata)

   var_names = strread(pinfo.var_names,'%s','delimiter',' ');

   nfigs = length(var_names);  % each variable gets its own figure
   iplot = 0;

   for ivar = 1:nfigs

      iplot = iplot + 1;
      figure(iplot); clf; 

      % Get the data from the netcdf files

      vname  = var_names{ivar};

      if ( strcmp(lower(vname),'ps') ==1 ) %  PS(time, copy, lat, lon)
         corner = [  1 NaN pinfo.latindex pinfo.lonindex];
         endpnt = [ -1 NaN pinfo.latindex pinfo.lonindex];
      else % U(time, copy, lat, lon, lev)
         corner = [  1 NaN pinfo.latindex pinfo.lonindex pinfo.levelindex];
         endpnt = [ -1 NaN pinfo.latindex pinfo.lonindex pinfo.levelindex];
      end

      ens_colors = get(gca,'ColorOrder');   % trying to cycle through colors
      ncolors = size(ens_colors,1) - 1;     % last one is black, already used.

      for i = 1:pinfo.copies

         imem = pinfo.copyindices(i);

         corner(2) = imem;
         endpnt(2) = imem;
         pr_series = getnc(pinfo.prior_file,     vname, corner, endpnt);
         po_series = getnc(pinfo.posterior_file, vname, corner, endpnt);

	 str1 = deblank(metadata(imem,:));

	 ens_member(1,:) = pr_series;
	 ens_member(2,:) = po_series;
	 b               = ens_member(:);

   	 hold on; 
	 memcolor = 1 + mod(i-1,ncolors); % cycles through colors [1,6]

         if (strcmp(str1,'ensemble mean') == 1)
	    h = plot(pinfo.xax, b,'k-','linewidth',2.0);
         else
	    h = plot(pinfo.xax, b,'linewidth',0.5,'Color',ens_colors(memcolor,:));
         end

	 [legh, objh, outh, outm] = legend;
	 nlines = length(outm);
         outm{nlines+1} = str1;
	 [legh, objh, outh, outm] = legend([outh; h],outm,0);
      end

      % Plot the true trajectory if it exists; the ens mean; annotate

      if ( exist(pinfo.truth_file) == 2 )
         corner(2) = truth.truth_index;
         endpnt(2) = truth.truth_index;
         true_trajectory = getnc(pinfo.truth_file, vname, corner, endpnt);
	 h = plot(prior.times, true_trajectory, 'k-','linewidth',1.0); hold on;
	 h = plot(pinfo.xax, a, 'k-','linewidth',2.0);
	 legend('truth','ensemble mean')
      end

      ylabel(sprintf('''%s''  (%.3f,%.3fE) level index %d', vname, ...
                     pinfo.latitude, pinfo.longitude, pinfo.levelindex ))
      xlabel(sprintf('model time (%d timesteps)', prior.num_times ))
      title(sprintf('%s Trajectories',pinfo.model), ...  
                     'interpreter','none','fontweight','bold')
      legend boxoff

   end
