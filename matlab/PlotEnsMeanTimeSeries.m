function PlotEnsMeanTimeSeries( pinfo )
%% PlotEnsMeanTimeSeries : plots time series of ensemble members, mean and truth, if the truth exists.
%
% PlotEnsMeanTimeSeries is intended to be called by 'plot_ens_mean_time_series'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE:    PlotEnsMeanTimeSeries( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
% var             name of netCDF variable of interest
% var_inds        indices of state variables of interest. Each variable gets
%                 plotted on its own axis.
%
%
% Example 1  (9variable model)
%-------------------------------------------------------------
% pinfo.truth_file  = 'True_State.nc';
% pinfo.diagn_file  = 'Prior_Diag.nc';
% pinfo.model       = '9var';
% pinfo.var         = 'state';
% pinfo.var_inds    = [ 4 5 6 ]; 

% PlotEnsMeanTimeSeries( pinfo )
%
% Example 2 (FMS BGrid model)
%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotEnsMeanTimeSeries( pinfo )

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( exist(pinfo.diagn_file,'file') ~= 2 ), error('%s does not exist.',pinfo.diagn_file); end

% Get the indices for the ensemble mean and spread.
% The metadata is queried to determine which "copy" is appropriate.
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% If the truth is known, great.
if ( exist(pinfo.truth_file,'file') == 2)
   have_truth  = 1;
   truth_index = get_copy_index(pinfo.truth_file, 'true state');
   vars        = CheckModelCompatibility(pinfo);
   pinfo       = CombineStructs(pinfo,vars);
else
   have_truth  = 0;
   diminfo = nc_getdiminfo(pinfo.diagn_file,'time');
   pinfo.diagn_time = [1 diminfo.Length];
   pinfo.model = nc_attget(pinfo.diagn_file, nc_global, 'model');
end

% Get some useful plotting arrays

num_times = pinfo.diagn_time(2) - pinfo.diagn_time(1) + 1;
times     = nc_varget(pinfo.diagn_file,'time', pinfo.diagn_time(1)-1, num_times); 

switch lower(pinfo.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3
            ivar = (i - 1)*3 + j;
            fprintf('plotting model %s Variable %d ...\n',pinfo.model,ivar)
            subplot(3, 1, j);

            % If the truth is known ...
            if ( have_truth )
               truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar, ...
                                      pinfo.truth_time(1), pinfo.truth_time(2));
               plot(times,truth,'b','LineWidth',1.0); hold on;
               legendstr = 'True State';
            end

            ens_mean = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar, ...
                                      pinfo.diagn_time(1), pinfo.diagn_time(2));
            plot(times,ens_mean,'r','LineWidth',1.0)

            if (exist('legendstr','var'))
               legend(legendstr,'Ensemble Mean',0)
            else
               legend('Ensemble Mean',0)
            end

            title(sprintf('%s Variable %d of %s',pinfo.model,ivar,pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',num_times))
            legend boxoff
            hold off;
         end
      end

   case {'lorenz_63','lorenz_84'}

      % Use one figure with three(usually) subplots
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,

            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);

            ens_mean    = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2));
            plot(times,ens_mean,'r');
            legend('Ensemble Mean',0);

            % If the truth is known ...
            if ( have_truth )
               truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar, ...
                                      pinfo.truth_time(1), pinfo.truth_time(2));
               hold on; plot(times,truth,'b'); hold off;
               legend('Ensemble Mean','True State',0);
            end

            title(sprintf('%s Variable %d of %s',pinfo.model,ivar,pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',num_times))
            legend boxoff
      end
      % as a bonus, plot the mean attractors.
      figure(2); clf
      ens  = get_state_copy(pinfo.diagn_file,pinfo.var, ens_mean_index, ...
                            pinfo.diagn_time(1), pinfo.diagn_time(2));
      plot3(ens(:,1), ens(:,2), ens(:,3), 'r')
      legend('Ensemble Mean',0)

      if (have_truth)
         ts= get_state_copy(pinfo.truth_file,pinfo.var, truth_index, ...
                            pinfo.truth_time(1), pinfo.truth_time(2));
         hold on; plot3(  ts(:,1),  ts(:,2),  ts(:,3), 'b'); hold off;
         legend('Ensemble Mean','True State',0)
      end

      title(sprintf('%s Attractors for %s', ...
              pinfo.model, pinfo.diagn_file), ...    
              'interpreter','none','fontweight','bold')
      xlabel('state variable 1')
      ylabel('state variable 2')
      zlabel('state variable 3')
      legend boxoff

   case {'lorenz_96', 'lorenz_96_2scale', 'lorenz_04', 'forced_lorenz_96', ...
         'ikeda', 'simple_advection'} 

      % Plot all variables in own subplot ... might get cluttered.
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);
            ens_mean = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2));
            plot(times,ens_mean,'r')
            legend('Ensemble Mean',0)

            % Get the truth for this variable
            if (have_truth)
               truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar, ...
                                      pinfo.truth_time(1), pinfo.truth_time(2));
               hold on; plot(times,truth,'b'); hold off;
               legend('Ensemble Mean','True State',0)
            end
            title(sprintf('%s Variable %d of %s',pinfo.model,ivar,pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',num_times))
            legend boxoff
      end

   case {'fms_bgrid','pe2lyr','mitgcm_ocean','wrf'}

      clf;

      timeunits = nc_attget(pinfo.fname, 'time',    'units');
      varunits  = nc_attget(pinfo.fname, pinfo.var, 'units');

      subplot(2,1,1)
         PlotLocator(pinfo)

      ens_mean    = GetCopy(pinfo.diagn_file, ens_mean_index,   pinfo , ...
                            pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

      subplot(2,1,2)
         plot(times,ens_mean,'r');
         legend('Ensemble Mean',0)

         s1 = sprintf('%s ''%s'' -- %s Ensemble Mean', ...
                            pinfo.model, pinfo.var, pinfo.diagn_file);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                    pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1,s2},'interpreter','none','fontweight','bold')

         if (have_truth)
            truth = GetCopy(pinfo.truth_file, truth_index,      pinfo , ...
                            pinfo.truth_time(1), pinfo.truth_time(2)) ;
            hold on; plot(times,truth,'b'); hold off;
            legend('Ensemble Mean','True State',0)
         end

         legend boxoff
         xlabel(sprintf('time (%s) %d timesteps',timeunits, num_times))
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

         timeunits = nc_attget(pinfo.diagn_file, 'time',    'units');
         varunits  = nc_attget(pinfo.diagn_file, pinfo.var, 'units');

         subplot(2,1,1)
            PlotLocator(pinfo);

         ens_mean = GetCopy(pinfo.diagn_file, ens_mean_index, pinfo, ...
                            pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
         subplot(2,1,2)
            plot(times,ens_mean,'r','LineWidth',2);
            legend('Ensemble Mean', 0)
            s1 = sprintf('%s model ''%s'' %s Ensemble Mean ', ...
                 pinfo.model, pinfo.var, pinfo.diagn_file);
            s2 = sprintf('level index %d lat %.2f lon %.2f', ...
                       pinfo.levelindex, pinfo.latitude, pinfo.longitude);

            if ( have_truth )
               truth    = GetCopy(pinfo.truth_file, truth_index, pinfo, ...
                                  pinfo.truth_time(1), pinfo.truth_time(2)) ;
               hold on; plot(times,truth,'b','LineWidth',2); hold off;
               legend('Ensemble Mean','True State',0);
               s1 = sprintf('%s model ''%s'' %s Truth and Ensemble Mean ', ...
                               pinfo.model, pinfo.var, pinfo.diagn_file);
            end

            %plot(times,ens_mean,'r','LineWidth',2); %      again - on top

            title({s1,s2},'interpreter','none','fontweight','bold')
            xlabel(sprintf('time (%s) %d timesteps',timeunits, num_times))
            ylabel(varunits)
            legend boxoff
      end

   otherwise
      error('model %s unknown.',pinfo.model)

end


%======================================================================
% Subfunctions
%======================================================================

function var = GetCopy(fname, copyindex, pinfo, tstartind, tendind)
% Gets a time-series of a single specified copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)
myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.levelindex = pinfo.levelindex;
myinfo.latindex   = pinfo.latindex;
myinfo.lonindex   = pinfo.lonindex;
[start, count]    = GetNCindices(myinfo,'diagn',pinfo.var);

varinfo = nc_getvarinfo(fname,pinfo.var);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tendind - tstartind + 1;
         break
      otherwise
   end
end
var = nc_varget(fname, pinfo.var, start, count);



function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axlims = axis;
   axlims = axlims + [-20 20 -20 20];
   grid on
   axis image
   axis(axlims)
   if (axlims(2) < 0)
       worldmap('hollow','dateline');
   else
       worldmap('hollow','greenwich');
   end

