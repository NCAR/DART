function PlotSawtooth( pinfo )
% PlotSawtooth: Plots time series of ensemble members, mean and truth.
%
%
% PlotSawtooth is intended to be called by 'plot_sawtooth'
% The only input argument is a structure with model-dependent
% components.
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
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html 

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

CheckModelCompatibility(pinfo.prior_file, pinfo.posterior_file)

% Get some information from the truth_file, if it exists.
if ( isfield(pinfo,'truth_file') ) 
   CheckModelCompatibility(pinfo.truth_file, pinfo.posterior_file)

   ft = netcdf(pinfo.truth_file);
   t.model      = ft.model(:);
   t.num_vars   = ncsize(ft('StateVariable')); % determine # of state variables
   t.num_copies = ncsize(ft('copy'));          % determine # of ensemble members
   t.num_times  = ncsize(ft('time'));          % determine # of output times
   close(ft);

   truth_index  = get_copy_index(pinfo.truth_file, 'true state' );
end

% Get some information from the prior_file 
fpr = netcdf(pinfo.prior_file);
pr.model      = fpr.model(:);
pr.num_vars   = ncsize(fpr('StateVariable')); % determine # of state variables
pr.num_copies = ncsize(fpr('copy'));          % determine # of ensemble members
pr.num_times  = ncsize(fpr('time'));          % determine # of output times
close(fpr);

% Get some information from the posterior_file 
fpo = netcdf(pinfo.posterior_file);
po.model      = fpo.model(:);
po.num_vars   = ncsize(fpo('StateVariable')); % determine # of state variables
po.num_copies = ncsize(fpo('copy'));          % determine # of ensemble members
po.num_times  = ncsize(fpo('time'));          % determine # of output times
po.timeunits  = fpo{'time'}.units(:);
po.calendar   = fpo{'time'}.calendar(:);
po.varunits   = fpo{pinfo.var}.units(:);
close(fpr);

% Get the indices for the true state and ensemble mean
% The metadata is queried to determine which "copy" is appropriate.
pr_ens_mean_index = get_copy_index(pinfo.prior_file, 'ensemble mean');
po_ens_mean_index = get_copy_index(pinfo.posterior_file, 'ensemble mean');

% Get some useful plotting arrays
times      = getnc(pinfo.prior_file,'time');
num_times  = length(times);
ens_mean   = zeros(2,num_times);   % trajectory of ensemble mean
ens_member = zeros(2,num_times);   % trajectory of ensemble member
x          = zeros(2,num_times);   % xaxis 
x(1,:)     = times;
x(2,:)     = times;
xax        = x(:);
metadata   = getnc(pinfo.prior_file,'CopyMetaData'); % get all the metadata

      nfigs = length(pinfo.var_inds);  % each variable gets its own figure
      iplot = 0;

      for ivar = pinfo.var_inds,

         % Get the data from the netcdf files

         po_ens_mean = get_var_series(pinfo.prior_file,    pinfo.var,po_ens_mean_index,ivar);
         pr_ens_mean = get_var_series(pinfo.posterior_file,pinfo.var,pr_ens_mean_index,ivar);


         % Now we paste them together in a clever way to show the effect of the assimilation

	 ens_mean(  1,:) = pr_ens_mean;
	 ens_mean(  2,:) = po_ens_mean;
	 a = ens_mean(:);       % unrolls them prior, post, prior, post ...

	 % Plot the true trajectory if it exists; the ens mean; annotate

         iplot = iplot + 1;
         figure(iplot); clf; 

	 if (isfield(pinfo,'truth_file'))
            truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
	    h = plot(times, truth,'k-','linewidth',1.0); hold on;
	    h = plot(xax, a, 'k-','linewidth',2.0);
	    legend('truth','ensemble mean')
         else 
	    h = plot(xax, a, 'k-','linewidth',2.0); hold on;
	    legend('ensemble mean')
         end

	 ylabel(sprintf('''%s'' index %d', pinfo.var, ivar ))
	 xlabel(sprintf('model time (%d timesteps)', num_times ))
	 title(sprintf('%s Trajectories',pinfo.model), ...  
                     'interpreter','none','fontweight','bold')

	 % Now check to see if we are overlaying any individual ensemble members.
	 % if pinfo.copyindices = [], nothing happens. 

	 ens_colors = get(gca,'ColorOrder');   % trying to cycle through colors
	 ncolors = size(ens_colors,1) - 1;     % last one is black, already used.

	 nmem = 0;
	 for imem = pinfo.copyindices
            nmem = nmem + 1;

            po_series   = get_var_series(pinfo.prior_file,     pinfo.var, imem, ivar);
            pr_series   = get_var_series(pinfo.posterior_file, pinfo.var, imem, ivar);

	    %str1 = sprintf('ensemble member %d',imem);
	    str1 = deblank(metadata(imem,:));

	    ens_member(1,:) = pr_series;
	    ens_member(2,:) = po_series;
	    b               = ens_member(:);     % unrolls them prior, post, prior, post ...

   	    hold on; 
	    memcolor = 1 + mod(nmem-1,ncolors); % cycles through colors [1,6]
	    h = plot(xax, b,'linewidth',0.5,'Color',ens_colors(memcolor,:));

	    [legh, objh, outh, outm] = legend;
	    nlines = length(outm);
            outm{nlines+1} = str1;
	    [legh, objh, outh, outm] = legend([outh; h],outm,0);
	 end
	 legend boxoff

      end
