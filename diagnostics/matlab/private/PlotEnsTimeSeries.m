function PlotEnsTimeSeries( pinfo )
%% PlotEnsTimeSeries: Plots time series of ensemble members, mean and truth.
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
% pinfo.truth_file  = 'true_state.nc';
% pinfo.diagn_file  = 'preassim.nc';
% pinfo.model       = '9var';
% pinfo.var         = 'state';
% pinfo.var_inds    = [ 4 5 6 ];
% PlotEnsTimeSeries( pinfo );
%
% Example 2 ( 9var model )
%%--------------------------------------------------------
% pinfo.truth_file  = 'true_state.nc';
% pinfo.diagn_file  = 'preassim.nc';
% pinfo.model       = 'fms_bgrid';
% pinfo.var         = 'u';
% pinfo.level       = 3;
% pinfo.latitude    = 23.5;
% pinfo.longitude   = 45.67;
% PlotEnsTimeSeries( pinfo )

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(pinfo.diagn_file,'file') ~= 2 ), error('%s does not exist.',pinfo.diagn_file); end

% If the truth is known, great.
if ( exist(pinfo.truth_file,'file') == 2)
   have_truth  = 1;
else
   have_truth  = 0;
   pinfo.diagn_time = [1 pinfo.time_series_length];
end

switch lower(pinfo.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3
            ivar = (i - 1)*3 + j;
            fprintf('plotting model %s Variable %d ...\n',pinfo.model,ivar)
            subplot(3, 1, j);

            ens_mean    = get_hyperslab('fname', pinfo.diagn_file, ...
                              'varname','state_mean', ...
                              'stateindex',ivar, ...
                              'tindex1',pinfo.diagn_time(1), ...
                              'tcount',pinfo.diagn_time(2)) ;

            ens_members = get_hyperslab('fname', pinfo.diagn_file, ...
                              'varname','state', ...
                              'stateindex',ivar, ...
                              'squeeze', 'true', ...
                              'tindex1', pinfo.diagn_time(1), ...
                              'tcount',pinfo.diagn_time(2)) ;

            h1 = plot(pinfo.time, ens_members, 'g'); hold on;
            h2 = plot(pinfo.time, ens_mean, 'r', 'LineWidth',2.0);
            h  = [h1(1), h2];
            legendstr = {sprintf('Ensemble Members (%d)',pinfo.num_ens_members), 'Ensemble Mean'};

            if (have_truth)
               truth    = get_hyperslab('fname',pinfo.truth_file, ...
                              'varname','state', ...
                              'stateindex',ivar, ...
                              'squeeze','true', ...
                              'tindex1', pinfo.truth_time(1),...
                              'tcount',pinfo.truth_time(2));

               h(3) = plot(pinfo.time, truth,'b','LineWidth',2.0);
               legendstr{3} = 'True State';
            end

            title(sprintf('%s Variable %d Ensemble Members of %s',...
                     pinfo.model, ivar, pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold');
            xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
            legend(h,legendstr)
            legend boxoff
            hold off;
         end
      end

   case {'lorenz_63','lorenz_84'}

      % Use one figure with three subplots
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,

            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);

            ens_mean    = get_hyperslab('fname',pinfo.diagn_file, ...
                              'varname','state_mean', ...
                              'stateindex',ivar, ...
                              'tindex1',pinfo.diagn_time(1), ...
                              'tcount',pinfo.diagn_time(2)) ;
            ens_members = get_hyperslab('fname',pinfo.diagn_file, ...
                              'varname','state', ...
                              'stateindex',ivar, ...
                              'squeeze','true', ...
                              'tindex1', pinfo.diagn_time(1), ...
                              'tcount',pinfo.diagn_time(2)) ;

            h1 = plot(pinfo.time,ens_members,'g'); hold on;
            h2 = plot(pinfo.time,   ens_mean,'r','LineWidth',2.0);
            h = [h1(1), h2];
            legendstr = {sprintf('Ensemble Members (%d)',pinfo.num_ens_members), 'Ensemble Mean'};

            if (have_truth)
               truth    = get_hyperslab('fname',pinfo.truth_file, ...
                              'varname','state', ...
                              'stateindex',ivar, ...
                              'squeeze','true',...
                              'tindex1', pinfo.truth_time(1),...
                               'tcount',pinfo.truth_time(2));
               h(3) = plot(pinfo.time, truth,'b','LineWidth',2.0);
               legendstr{3} = 'True State';
            end

            title(sprintf('%s Variable %d Ensemble Members of %s', ...
                  pinfo.model, ivar, pinfo.diagn_file), ...
                  'interpreter','none','fontweight','bold');
            xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
            legend(h,legendstr)
            legend boxoff
            hold off;
      end

      % as a bonus, plot the mean attractors.
      figure(2); clf

      if (have_truth)
         ts   = get_hyperslab('fname',pinfo.truth_file, ...
                              'varname','state', ...
                              'tindex1',pinfo.diagn_time(1), ...
                              'tcount',pinfo.diagn_time(2));
         plot3( ts(1,:), ts(2,:), ts(3,:), 'b'); hold on;
         legendstr = 'True State';
      end

      ens  = get_hyperslab('fname',pinfo.diagn_file, ...
                           'varname','state_mean', ...
                           'tindex1',pinfo.diagn_time(1), ...
                           'tcount',pinfo.diagn_time(2));
      plot3(ens(1,:), ens(2,:), ens(3,:), 'r');
      title(sprintf('%s Attractors for %s', pinfo.model, pinfo.diagn_file), ...
                 'interpreter','none','fontweight','bold');

      if (exist('legendstr','var'))
         legend(legendstr,'Ensemble Mean','Location','NorthEast');
      else
         legend(          'Ensemble Mean','Location','NorthEast');
      end

      xlabel('state variable 1');
      ylabel('state variable 2');
      zlabel('state variable 3');
      legend boxoff

   case {'lorenz_96', 'lorenz_96_2scale', 'forced_lorenz_96', 'lorenz_04', ...
         'ikeda', 'simple_advection', 'lorenz_96_tracer_advection', 'null'}

      % Use one figure with subplots
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);

            varname = sprintf('%s_mean',pinfo.var);
            ens_mean    = get_hyperslab('fname',pinfo.diagn_file, ...
                              'varname',varname, ...
                              'stateindex',ivar, ...
                              'tindex1',pinfo.diagn_time(1), ...
                              'tcount',pinfo.diagn_time(2)) ;

            ens_members = get_hyperslab('fname',pinfo.diagn_file, ...
                              'varname',pinfo.var, ...
                              'stateindex',ivar, ...
                              'squeeze', 'true', ...
                              'tindex1', pinfo.diagn_time(1), ...
                              'tcount',pinfo.diagn_time(2)) ;

            h1 = plot(pinfo.time,ens_members,'g'); hold on;
            h2 = plot(pinfo.time,   ens_mean,'r','LineWidth',2.0);
            h = [h1(1), h2];
            legendstr = {sprintf('Ensemble Members (%d)',pinfo.num_ens_members), 'Ensemble Mean'};

            if (have_truth)
               truth    = get_hyperslab('fname',pinfo.truth_file, ...
                              'varname',pinfo.var, ...
                              'stateindex',ivar, ...
                              'squeeze', 'true', ...
                              'tindex1', pinfo.truth_time(1), ...
                              'tcount',pinfo.truth_time(2));
               h(3) = plot(pinfo.time, truth,'b','LineWidth',2.0);
               legendstr{3} = 'True State';
            end

            title(sprintf('%s %s varnum %d Ensemble Members of %s',...
                     pinfo.model, pinfo.var, ivar, pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
            xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
            legend(h,legendstr)
            legend boxoff
            hold off;
      end

   case {'fms_bgrid','pe2lyr','mitgcm_ocean','wrf','cam','sqg','pop'}

      clf;

      varunits  = ncreadatt(pinfo.fname, pinfo.var, 'units');

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)

         varname = sprintf('%s_mean',pinfo.var);
         ens_mean = get_hyperslab('fname',pinfo.diagn_file, ...
                   'varname',varname, ...
                   'permute','T', ...
                   'levelindex',pinfo.levelindex, ...
                   'latindex',pinfo.latindex, ...
                   'lonindex',pinfo.lonindex, ...
                   'tindex1',pinfo.diagn_time(1), ...
                   'tcount',pinfo.diagn_time(2));

         ens_members = get_hyperslab('fname',pinfo.diagn_file, ...
                   'varname',pinfo.var, ...
                   'permute','T', ...
                   'levelindex',pinfo.levelindex, ...
                   'latindex',pinfo.latindex, ...
                   'lonindex',pinfo.lonindex, ...
                   'tindex1',pinfo.diagn_time(1), ...
                   'tcount',pinfo.diagn_time(2));

         hmems = plot(pinfo.time, ens_members,'g-'); hold on;
         hmean = plot(pinfo.time, ens_mean,   'r-','LineWidth',2);
         legendstr = {sprintf('Ensemble Members (%d)',pinfo.num_ens_members),...
                      'Ensemble Mean'};
         h = [hmems(1), hmean];

         if ( have_truth )
            truth = get_hyperslab('fname',pinfo.truth_file, ...
                   'varname',pinfo.var, ...
                   'permute','T', ...
                   'levelindex',pinfo.levelindex, ...
                   'latindex',pinfo.latindex, ...
                   'lonindex',pinfo.lonindex, ...
                   'tindex1',pinfo.truth_time(1), ...
                   'tcount',pinfo.truth_time(2));

            h(3) = plot(pinfo.time, truth, 'b-','LineWidth',2);
            legendstr{3} = 'True State';
         end
         s1 = sprintf('%s model ''%s'' %s Ensemble Members ', ...
                    pinfo.model, pinfo.var, pinfo.diagn_file);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                    pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1,s2},'interpreter','none','fontweight','bold');
         xdates(pinfo.time)
         ylabel(varunits);
         legend(h,legendstr);
         legend boxoff
         hold off;

    case {'mpas_atm'}

      clf;

      varunits  = ncreadatt(pinfo.fname, pinfo.var, 'units');

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)

         varname = sprintf('%s_mean',pinfo.var);
         ens_mean = get_hyperslab('fname',pinfo.diagn_file, ...
                   'varname',varname, ...
                   'permute','T', ...
                   'levelindex',pinfo.levelindex, ...
                   'copyindex',ens_mean_index, ...
                   'cellindex',pinfo.cellindex, ...
                   'tindex1',pinfo.diagn_time(1), ...
                   'tcount',pinfo.diagn_time(2));

         ens_members = get_hyperslab('fname',pinfo.diagn_file, ...
                   'varname',pinfo.var, ...
                   'permute','T', ...
                   'levelindex',pinfo.levelindex, ...
                   'cellindex',pinfo.cellindex, ...
                   'tindex1',pinfo.diagn_time(1), ...
                   'tcount',pinfo.diagn_time(2));

         hmems = plot(pinfo.time, ens_members,'g-'); hold on;
         hmean = plot(pinfo.time, ens_mean,   'r-','LineWidth',2);
         legendstr = {sprintf('Ensemble Members (%d)',pinfo.num_ens_members),...
                      'Ensemble Mean'};
         h = [hmems(1) hmean];

         if ( have_truth )
            truth  = get_hyperslab('fname',pinfo.truth_file, ...
                   'varname',pinfo.var, ...
                   'permute','T', ...
                   'levelindex',pinfo.levelindex, ...
                   'cellindex',pinfo.cellindex, ...
                   'tindex1',pinfo.truth_time(1), ...
                   'tcount',pinfo.truth_time(2));

            h(3) = plot(pinfo.time, truth, 'b-','LineWidth',2);
            legendstr{3} = 'True State';
         end

         s1 = sprintf('%s model ''%s'' %s Ensemble Members ', ...
                    pinfo.model, pinfo.var, pinfo.diagn_file);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                    pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1,s2},'interpreter','none','fontweight','bold');
         xdates(pinfo.time)
         ylabel(varunits);
         legend(h,legendstr);
         legend boxoff

   otherwise

      error('model %s unknown.',pinfo.model)

end


%======================================================================
% Subfunctions
%======================================================================



function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axlims = axis;
   axlims = axlims + [-20 20 -20 20];
   grid on
   axis image
   axis(axlims)
   if (axlims(2) < 0)
       continents('hollow','dateline');
   else
       continents('hollow','greenwich');
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
