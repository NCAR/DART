function PlotEnsTimeSeries( pinfo )
% PlotEnsTimeSeries: Plots time series of ensemble members, mean and truth.
%
%
% PlotEnsTimeSeries is intended to be called by 'plot_ens_time_series'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotEnsTimeSeries( pinfo );
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
% var             name of netCDF variable of interest
% var_inds        indices of variables of interest
%
% Example 1 ( 9var model )
%%--------------------------------------------------------
% pinfo.truth_file  = 'True_State.nc';
% pinfo.diagn_file  = 'Posterior_Diag.nc';
% pinfo.var         = 'state';
% pinfo.var_inds    = [ 4 5 6 ]; 
% PlotEnsTimeSeries( pinfo );
%
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotEnsTimeSeries( pinfo )

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


if ( exist(pinfo.truth_file) == 2)
   CheckModelCompatibility(pinfo.truth_file, pinfo.diagn_file)
   truth_index = get_copy_index(pinfo.truth_file, 'true state' );
   have_truth  = 1;
else
   have_truth  = 0;
end

% Get some information from the diagn_file 
fd = netcdf(pinfo.diagn_file);
d.model      = fd.model(:);
d.num_copies = length(fd('copy')); % determine # of ensemble members
d.num_times  = length(fd('time')); % determine # of output times
close(fd);

% Get the indices for the ensemble mean
% The metadata is queried to determine which "copy" is appropriate.
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');

% Get some useful plotting arrays
times = getnc(pinfo.diagn_file,'time');

switch lower(d.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3

            ivar = (i - 1)*3 + j;
            ens_mean    = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar);
            ens_members = get_ens_series(pinfo.diagn_file, pinfo.var, ivar);  % all members
            if (have_truth) 
               truth    = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
            end

            subplot(3, 1, j);
               % This is a bit wasteful, but we plot everything once to define
               % the color order for the legend and then again for visibility
               plot(times,      truth,'b','LineWidth',2.0); hold on;
               plot(times,   ens_mean,'r','LineWidth',2.0); hold on;
               plot(times,ens_members,'g');
               legend('True State','Ensemble Mean', ...
                      sprintf('Ensemble Members (%d)',d.num_copies-2),0)

               legend boxoff
               plot(times,   truth,'b','LineWidth',2.0); % plot again - on top
               plot(times,ens_mean,'r','LineWidth',2.0); %      again - on top
               title(sprintf('%s Variable %d Ensemble Members of %s',...
                     d.model, ivar, pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
               xlabel(sprintf('model time (%d timesteps)',d.num_times))
         end
      end

   case {'lorenz_63','lorenz_84'}

      % Use one figure with three subplots 
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,

            truth       = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
            ens_mean    = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar);
            ens_members = get_ens_series(pinfo.diagn_file, pinfo.var, ivar);  % all members

            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);
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
                     d.model, ivar, pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
               xlabel(sprintf('model time (%d timesteps)',d.num_times))
      end
      % as a bonus, plot the mean attractors.
      figure(2); clf
      ts   = get_state_copy(pinfo.truth_file,pinfo.var, truth_index);
      ens  = get_state_copy(pinfo.diagn_file,pinfo.var, ens_mean_index);
      plot3(  ts(:,1),  ts(:,2),  ts(:,3), 'b', ...
             ens(:,1), ens(:,2), ens(:,3), 'r')
      title(sprintf('%s Attractors for %s and %s', ...
                 d.model, pinfo.truth_file, pinfo.diagn_file), ...    
                 'interpreter','none','fontweight','bold')
      legend('True State','Ensemble Mean',0)
      legend boxoff
      xlabel('state variable 1')
      ylabel('state variable 2')
      zlabel('state variable 3')

   case {'lorenz_96', 'lorenz_96_2scale', 'forced_lorenz_96', 'lorenz_04'}

      % Use one figure with subplots 
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,

            truth       = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
            ens_mean    = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar);
            ens_members = get_ens_series(pinfo.diagn_file, pinfo.var, ivar);  % all members

            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);
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
               title(sprintf('%s %s varnum %d Ensemble Members of %s',...
                     d.model, pinfo.var, ivar, pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
               xlabel(sprintf('model time (%d timesteps)',d.num_times))
      end

   case 'fms_bgrid'

      clf;

      ft = netcdf(pinfo.truth_file);
      d.model      = fd.model(:);
      timeunits    = ft{'time'}.units(:);
      varunits     = ft{pinfo.var}.units(:);
      close(ft);

      truth       = GetCopy(pinfo.truth_file, truth_index,      pinfo );
      ens_mean    = GetCopy(pinfo.diagn_file, ens_mean_index,   pinfo );
      ens_members = GetEns( pinfo.diagn_file,                   pinfo );

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)

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

         s1 = sprintf('%s model ''%s'' Truth and %s Ensemble Members ', ...
                            d.model, pinfo.var, pinfo.diagn_file);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                    pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1,s2},'interpreter','none','fontweight','bold')

         xlabel(sprintf('time (%s) %d timesteps',timeunits, d.num_times))
         ylabel(varunits)

   case 'cam'

      clf;

      var_names = strread(pinfo.var_names,'%s','delimiter',' ');
      nfigs = length(var_names);  % each variable gets its own figure
      iplot = 0;

      for ivar = 1:nfigs

         iplot = iplot + 1;
         figure(iplot); clf;

         pinfo.var  = var_names{ivar};

         fd         = netcdf(pinfo.diagn_file);
         model      = fd.model(:);
         timeunits  = fd{'time'}.units(:);
         varunits   = fd{pinfo.var}.units(:);
         close(fd);

         if ( have_truth )
            truth    = GetCamCopy(pinfo.truth_file, truth_index,      pinfo );
         end

         ens_mean    = GetCamCopy(pinfo.diagn_file, ens_mean_index,   pinfo );
         ens_members = GetCamEns( pinfo.diagn_file,                   pinfo );
         nmembers    = size(ens_members,2);

         subplot(2,1,1)
            PlotLocator(pinfo);

         subplot(2,1,2)
            % This is a bit wasteful, but we plot everything once to define
            % the color order for the legend and then again for visibility
            if (have_truth)
               plot(times,      truth,'b','LineWidth',2); hold on;
               plot(times,   ens_mean,'r','LineWidth',2);
               plot(times,ens_members,'g');
               legend('True State','Ensemble Mean', ...
                            sprintf('Ensemble Members (%d)',nmembers),0);
               plot(times,   truth,'b','LineWidth',2); % plot again - on top
               s1 = sprintf('%s model ''%s'' Truth and %s Ensemble Members ', ...
                               model, pinfo.var, pinfo.diagn_file);
            else
               plot(times,   ens_mean,'r','LineWidth',2); hold on;
               plot(times,ens_members,'g');
               legend('Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers),0)
               s1 = sprintf('%s model ''%s'' %s Ensemble Members ', ...
                               model, pinfo.var, pinfo.diagn_file);
            end
            legend boxoff
            plot(times,ens_mean,'r','LineWidth',2); %      again - on top

            s2 = sprintf('level index %d lat %.2f lon %.2f', ...
                       pinfo.levelindex, pinfo.latitude, pinfo.longitude);
            title({s1,s2},'interpreter','none','fontweight','bold')

            xlabel(sprintf('time (%s) %d timesteps',timeunits, d.num_times))
            ylabel(varunits)
      end

   otherwise

      error(sprintf('model %s unknown.',d.model))

end


%======================================================================
% Subfunctions
%======================================================================

function var = GetCopy(fname, copyindex, pinfo)
% Gets a time-series of a single specified copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)
if strcmp(lower(pinfo.var),'ps')
   corner = [ 1 copyindex                  pinfo.latindex pinfo.lonindex];
   endpnt = [-1 copyindex                  pinfo.latindex pinfo.lonindex];
else
   corner = [ 1 copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
   endpnt = [-1 copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
end
var = getnc(fname, pinfo.var, corner, endpnt);



function var = GetCamCopy(fname, copyindex, pinfo)
% Gets a time-series of a single specified copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)
if strcmp(lower(pinfo.var),'ps')
   corner = [ 1 copyindex pinfo.latindex pinfo.lonindex];
   endpnt = [-1 copyindex pinfo.latindex pinfo.lonindex];
else
   corner = [ 1 copyindex pinfo.latindex pinfo.lonindex pinfo.levelindex];
   endpnt = [-1 copyindex pinfo.latindex pinfo.lonindex pinfo.levelindex];
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
if strcmp(lower(pinfo.var),'ps')
   corner = [-1 -1                  pinfo.latindex pinfo.lonindex];
   endpnt = [-1 -1                  pinfo.latindex pinfo.lonindex];
else
   corner = [-1 -1 pinfo.levelindex pinfo.latindex pinfo.lonindex];
   endpnt = [-1 -1 pinfo.levelindex pinfo.latindex pinfo.lonindex];
end
bob = getnc(fname, pinfo.var, corner, endpnt); % 'bob' is only 2D 
var = bob(:,copyindices);



function var = GetCamEns(fname, pinfo)
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
if strcmp(lower(pinfo.var),'ps')
   corner = [-1 -1 pinfo.latindex pinfo.lonindex];
   endpnt = [-1 -1 pinfo.latindex pinfo.lonindex];
else
   corner = [-1 -1 pinfo.latindex pinfo.lonindex pinfo.levelindex];
   endpnt = [-1 -1 pinfo.latindex pinfo.lonindex pinfo.levelindex];
end

bob = getnc(fname, pinfo.var, corner, endpnt); % 'bob' is only 2D 
var = bob(:,copyindices);



function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pg','MarkerSize',12,'MarkerFaceColor','g');
   axis([0 360 -90 90])
   worldmap;
   axis image
   grid on

