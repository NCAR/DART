function PlotEnsMeanTimeSeries( pinfo )
% PlotEnsMeanTimeSeries : plots time series of ensemble members, mean and truth 
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
%%-------------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'state';
% pinfo.var_inds   = [ 4 5 6 ];
% PlotEnsMeanTimeSeries( pinfo )
%
% Example 2 (FMS BGrid model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotEnsMeanTimeSeries( pinfo )

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

CheckModelCompatibility(pinfo.truth_file, pinfo.diagn_file)

% Get some information from the truth_file 
ft = netcdf(pinfo.truth_file);
t.model      = ft.model(:);
tvar_atts    = dim(ft{pinfo.var});      % cell array of dimensions for the var
t.num_times  = length(tvar_atts{1});    % determine # of output times
t.num_copies = length(tvar_atts{2});    % # of ensemble members
t.num_vars   = length(tvar_atts{3});    % dimension of desired variable
close(ft);

% Get some information from the diagn_file 
fd = netcdf(pinfo.diagn_file);
d.model      = fd.model(:);
dvar_atts    = dim(fd{pinfo.var});      % cell array of dimensions for the var
d.num_times  = length(dvar_atts{1});    % determine # of output times
d.num_copies = length(dvar_atts{2});    % # of ensemble members
d.num_vars   = length(dvar_atts{3});    % dimension of desired variable
close(fd);

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Get the appropriate copies
% (TJH) This function gets them "on the fly" in the plotting.
%truth      = get_state_copy(pinfo.truth_file, pinfo.var, truth_index);
%ens_mean   = get_state_copy(pinfo.diagn_file, pinfo.var, ens_mean_index );
%ens_spread = get_state_copy(pinfo.diagn_file, pinfo.var, ens_spread_index );

% Get some useful plotting arrays
times = getnc(pinfo.truth_file,'time');

switch lower(t.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3,
         figure(i); clf
         for j = 1:3,
            ivar = (i - 1)*3 + j;
            disp(sprintf('plotting model %s Variable %d ...',t.model,ivar))
            % Get the truth for this variable
            truth    = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
            ens_mean = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar);
            subplot(3, 1, j);
            plot(times,truth, 'b',times,ens_mean,'r')
            title(sprintf('%s Variable %d of %s',t.model,ivar,pinfo.diagn_file), ...
                  'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',t.num_times))
            legend('True State','Ensemble Mean',0)
            legend boxoff
         end
      end

   case {'lorenz_63','lorenz_84'}

      % Use one figure with three(usually) subplots
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            % Get the truth for this variable
            truth    = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
            ens_mean = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar);
            subplot(length(pinfo.var_inds), 1, iplot);
            plot(times,truth, 'b',times,ens_mean,'r')
            title(sprintf('%s Variable %d of %s',t.model,ivar,pinfo.diagn_file), ...
                  'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',t.num_times))
            legend('True State','Ensemble Mean',0)
            legend boxoff
      end
      % as a bonus, plot the mean attractors.
      figure(2); clf
      ts   = get_state_copy(pinfo.diagn_file, pinfo.var, truth_index);
      ens  = get_state_copy(pinfo.diagn_file, pinfo.var, ens_mean_index);
      plot3(  ts(:,1),  ts(:,2),  ts(:,3), 'b', ...
             ens(:,1), ens(:,2), ens(:,3), 'r')
      title(sprintf('%s Attractors for %s and %s', ...
                 t.model, pinfo.truth_file, pinfo.diagn_file), ...
                 'interpreter','none','fontweight','bold')
      legend('True State','Ensemble Mean',0)
      legend boxoff
      xlabel('state variable 1')
      ylabel('state variable 2')
      zlabel('state variable 3')

   case {'lorenz_96', 'lorenz_96_2scale', 'lorenz_04', 'forced_lorenz_96'} 
      
      % Plot all variables in own subplot ... might get cluttered.
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            % Get the truth for this variable
            truth    = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
            ens_mean = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar);
            subplot(length(pinfo.var_inds), 1, iplot);
            plot(times,truth, 'b',times,ens_mean,'r')
            title(sprintf('%s Variable %d of %s',t.model,ivar,pinfo.diagn_file), ...
                  'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',t.num_times))
            legend('True State','Ensemble Mean',0)
            legend boxoff
      end

   case 'fms_bgrid'

      clf;

      % Get some plotting information from the truth_file 
      ft = netcdf(pinfo.truth_file);
      timeunits = ft{'time'}.units(:);
      varunits  = ft{pinfo.var}.units(:);
      close(ft);

      truth      = GetCopy(pinfo.truth_file, truth_index,      pinfo );
      ens_mean   = GetCopy(pinfo.diagn_file, ens_mean_index,   pinfo );

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
         plot(times,truth, 'b', times,ens_mean, 'r');

         s1 = sprintf('%s ''%s'' -- Truth and %s Ensemble Mean', ...
                            t.model, pinfo.var, pinfo.diagn_file);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                    pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1,s2},'interpreter','none','fontweight','bold')

         legend('True State','Ensemble Mean',0)
         legend boxoff
         xlabel(sprintf('time (%s) %d timesteps',timeunits, t.num_times))
         ylabel(varunits)

   otherwise
      error(sprintf('model %s unknown.',t.model))

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


function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pg','MarkerSize',12,'MarkerFaceColor','g');
   axis([0 360 -90 90])
   worldmap
   axis image
   grid on


