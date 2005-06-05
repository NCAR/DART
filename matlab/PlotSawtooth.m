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

CheckModelCompatibility(pinfo.truth_file, pinfo.prior_file)
CheckModelCompatibility(pinfo.truth_file, pinfo.posterior_file)

% Get some information from the truth_file 
ft = netcdf(pinfo.truth_file);
t.model      = ft.model(:);
t.num_vars   = ncsize(ft('StateVariable')); % determine # of state variables
t.num_copies = ncsize(ft('copy')); % determine # of ensemble members
t.num_times  = ncsize(ft('time')); % determine # of output times
close(ft);

% Get some information from the prior_file 
fpr = netcdf(pinfo.prior_file);
pr.model      = fpr.model(:);
pr.num_vars   = ncsize(fpr('StateVariable')); % determine # of state variables
pr.num_copies = ncsize(fpr('copy')); % determine # of ensemble members
pr.num_times  = ncsize(fpr('time')); % determine # of output times
close(fpr);

% Get some information from the posterior_file 
fpo = netcdf(pinfo.posterior_file);
po.model      = fpo.model(:);
po.num_vars   = ncsize(fpo('StateVariable')); % determine # of state variables
po.num_copies = ncsize(fpo('copy')); % determine # of ensemble members
po.num_times  = ncsize(fpo('time')); % determine # of output times
close(fpr);

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index       = get_copy_index(pinfo.truth_file, 'true state' );
pr_ens_mean_index = get_copy_index(pinfo.prior_file, 'ensemble mean');
pr_ens_sprd_index = get_copy_index(pinfo.prior_file, 'ensemble spread');
po_ens_mean_index = get_copy_index(pinfo.posterior_file, 'ensemble mean');
po_ens_sprd_index = get_copy_index(pinfo.posterior_file, 'ensemble spread');

% Get some useful plotting arrays
times = getnc(pinfo.truth_file,'time');

switch lower(t.model)

   case 'fms_bgrid'

      clf;

      ft = netcdf(pinfo.truth_file);
      t.model      = ft.model(:);
      timeunits    = ft{'time'}.units(:);
      varunits     = ft{pinfo.var}.units(:);
      close(ft);

      truth       = GetCopy(pinfo.truth_file, truth_index,      pinfo );
      ens_mean    = GetCopy(pinfo.prior_file, ens_mean_index,   pinfo );
      ens_members = GetEns( pinfo.prior_file,                   pinfo );

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)

         % This is a bit wasteful, but we plot everything once to define
         % the color order for the legend and then again for visibility
         plot(times,      truth,'b','LineWidth',2); hold on;
         plot(times,   ens_mean,'r','LineWidth',2);
         plot(times,ens_members,'g');
         legend('True State','Ensemble Mean', ...
                      sprintf('Ensemble Members (%d)',pr.num_copies-2),0)
         legend boxoff
         plot(times,   truth,'b','LineWidth',2); % plot again - on top
         plot(times,ens_mean,'r','LineWidth',2); %      again - on top

         s1 = sprintf('%s model ''%s'' Truth and %s Ensemble Members ', ...
                            t.model, pinfo.var, pinfo.prior_file);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                    pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1,s2},'interpreter','none','fontweight','bold')

         xlabel(sprintf('time (%s) %d timesteps',timeunits, t.num_times))
         ylabel(varunits)

    otherwise

      nfigs = length(pinfo.var_inds);  % each variable gets its own figure
      iplot = 0;

      ens_mean   = zeros(2,length(times));   % trajectory of ensemble mean
      ens_member = zeros(2,length(times));   % trajectory of ensemble member
      x          = zeros(2,length(times));   % xaxis 
      x(1,:) = times;
      x(2,:) = times;
      xax    = x(:);

      for ivar = pinfo.var_inds,

         truth       = get_var_series(pinfo.truth_file,     pinfo.var, truth_index,       ivar);
         po_ens_mean = get_var_series(pinfo.prior_file,     pinfo.var, po_ens_mean_index, ivar);
         po_ens_sprd = get_var_series(pinfo.prior_file,     pinfo.var, po_ens_sprd_index, ivar);
         pr_ens_mean = get_var_series(pinfo.posterior_file, pinfo.var, pr_ens_mean_index, ivar);
         pr_ens_sprd = get_var_series(pinfo.posterior_file, pinfo.var, pr_ens_sprd_index, ivar);

	 str1 = sprintf('ensemble member %d',pinfo.copyindices(1));

         po_copyindx = get_copy_index(pinfo.prior_file,     str1);
         po_series   = get_var_series(pinfo.prior_file,     pinfo.var, po_copyindx, ivar);

         pr_copyindx = get_copy_index(pinfo.posterior_file, str1);
         pr_series   = get_var_series(pinfo.posterior_file, pinfo.var, pr_copyindx, ivar);

         % Now we paste them together in a clever way to show the effect of the assimilation

	 ens_mean(  1,:) = pr_ens_mean;
	 ens_mean(  2,:) = po_ens_mean;
	 ens_member(1,:) = pr_series;
	 ens_member(2,:) = po_series;

	 a = ens_mean(:);       % unrolls them prior, post, prior, post ...
	 b = ens_member(:);     % unrolls them prior, post, prior, post ...

	 plot(times, truth,'k', xax, a, 'r', xax, b, 'c')
	 legend('truth','ensemble mean',sprintf('ensemble member %2d',pinfo.copyindices(1)))
	 legend boxoff

	 break

         figure(1); clf; 

            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);
               % This is a bit wasteful, but we plot everything once to define
               % the color order for the legend and then again for visibility
               plot(times,      truth,'b','LineWidth',1.0); hold on;
               plot(times,   ens_mean,'r','LineWidth',1.0);
               plot(times,ens_members,'g');
               legend('True State','Ensemble Mean', ...
                      sprintf('Ensemble Members (%d)',pr.num_copies-2),0)
               legend boxoff
               plot(times,   truth,'b','LineWidth',1.0); % plot again - on top
               plot(times,ens_mean,'r','LineWidth',1.0); %      again - on top
               title(sprintf('%s Variable %d Ensemble Members of %s',...
                     t.model, ivar, pinfo.prior_file), ...
                     'interpreter','none','fontweight','bold')
               xlabel(sprintf('model time (%d timesteps)',t.num_times))
      end

end


%======================================================================
% Subfunctions
%======================================================================

function var = GetCopy(fname, copyindex, pinfo)
% Gets a time-series of a single specified copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)
if strcmp(pinfo.var,'ps')
   corner = [-1 copyindex                  pinfo.latindex pinfo.lonindex];
   endpnt = [-1 copyindex                  pinfo.latindex pinfo.lonindex];
else
   corner = [-1 copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
   endpnt = [-1 copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
end
var = getnc(fname, pinfo.var, corner, endpnt);



function var = GetEns(fname, pinfo)
% Gets a time-series of all copies of a prognostic variable 
% at a particular 3D location (level, lat, lon).
% Determining just the ensemble members (and not mean, spread ...)
% is the hard part.

% find which are actual ensemble members
metadata    = getnc(fname,'CopyMetaData');           % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's

if ( isempty(copyindices) )
   disp(sprintf('%s has no valid ensemble members',fname))
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   disp(sprintf('%s claims to have %d copies',fname, num_copies))
   error('netcdf file has no ensemble members.')
end
ens_num     = length(copyindices);

% Get all ensemble members, just return desired ones.
if strcmp(pinfo.var,'ps')
   corner = [-1 -1                  pinfo.latindex pinfo.lonindex];
   endpnt = [-1 -1                  pinfo.latindex pinfo.lonindex];
else
   corner = [-1 -1 pinfo.levelindex pinfo.latindex pinfo.lonindex];
   endpnt = [-1 -1 pinfo.levelindex pinfo.latindex pinfo.lonindex];
end
bob = getnc(fname, pinfo.var, corner, endpnt); % 'bob' is only 2D 
var = bob(:,copyindices);



function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pg','MarkerSize',12,'MarkerFaceColor','g');
   axis([0 360 -90 90])
   worldmap
   axis image
   grid on



