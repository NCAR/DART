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
% pinfo.model       = '9var';
% pinfo.var         = 'state';
% pinfo.var_inds    = [ 4 5 6 ]; 
% PlotEnsTimeSeries( pinfo );
%
% Example 2 ( 9var model )
%%--------------------------------------------------------
% pinfo.truth_file  = 'True_State.nc';
% pinfo.diagn_file  = 'Prior_Diag.nc';
% pinfo.model       = 'fms_bgrid';
% pinfo.var         = 'u';
% pinfo.level       = 3;
% pinfo.latitude    = 23.5;
% pinfo.longitude   = 45.67;
% PlotEnsTimeSeries( pinfo )

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( exist(pinfo.diagn_file,'file') ~= 2 ), error('%s does not exist.',pinfo.diagn_file); end

% Get the indices for the ensemble mean.
% The metadata is queried to determine which "copy" is appropriate.
ens_mean_index = get_copy_index(pinfo.diagn_file, 'ensemble mean');

% If the truth is known, great.
if ( exist(pinfo.truth_file,'file') == 2)
   have_truth  = 1;
   truth_index = get_copy_index(pinfo.truth_file, 'true state' );
   vars        = CheckModelCompatibility(pinfo);
   pinfo       = CombineStructs(pinfo,vars);
   pinfo.truth_times = nc_varget(pinfo.truth_file, 'time');
else
   have_truth  = 0;
   diminfo = nc_getdiminfo(pinfo.diagn_file,'time');
   pinfo.diagn_time = [1 diminfo.Length];
   pinfo.model = nc_attget(pinfo.diagn_file, nc_global, 'model');
end

% Get some useful plotting arrays
num_times = pinfo.diagn_time(2) - pinfo.diagn_time(1) + 1;
times     = nc_varget(pinfo.diagn_file,'time', pinfo.diagn_time(1)-1, num_times);

% Get some information from the diagn_file 
d.num_copies = dim_length(pinfo.diagn_file,'copy');

switch lower(pinfo.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3
            ivar = (i - 1)*3 + j;
            fprintf('plotting model %s Variable %d ...\n',pinfo.model,ivar)
            subplot(3, 1, j);

            if (have_truth) 
               truth    = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
               plot(pinfo.truth_times, truth,'b','LineWidth',2.0); hold on;
               legendstr = 'True State';
            end

            ens_mean    = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
            ens_members = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
            nmembers    = size(ens_members,2);

            plot(times,   ens_mean,'r','LineWidth',2.0); hold on;
            plot(times,ens_members,'g');
            if (exist('legendstr')) 
                legend(legendstr, 'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
               plot(pinfo.truth_times, truth,'b','LineWidth',2); % again, to put 'on top'
            else
               legend(           'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
            end
            plot(times,   ens_mean,'r','LineWidth',2.0);

            title(sprintf('%s Variable %d Ensemble Members of %s',...
                     pinfo.model, ivar, pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold');
            xlabel(sprintf('model time (%d timesteps)',num_times));
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

            if (have_truth)
               truth    = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
               plot(pinfo.truth_times, truth,'b','LineWidth',2.0); hold on;
               legendstr = 'True State';
            end

            ens_mean    = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2));
            ens_members = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2));
            nmembers    = size(ens_members,2);

            plot(times,   ens_mean,'r','LineWidth',2.0); hold on;
            plot(times,ens_members,'g');
            if (exist('legendstr'))
               legend(legendstr, 'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
               plot(pinfo.truth_times,   truth,'b','LineWidth',2); % again, to put 'on top'
            else
               legend(           'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
            end
            plot(times,   ens_mean,'r','LineWidth',2.0); % again, to put 'on top'

            title(sprintf('%s Variable %d Ensemble Members of %s', ...
                  pinfo.model, ivar, pinfo.diagn_file), ...
                  'interpreter','none','fontweight','bold');
            xlabel(sprintf('model time (%d timesteps)',num_times));
            legend boxoff
            hold off;
      end

      % as a bonus, plot the mean attractors.
      figure(2); clf

      if (have_truth)
         ts   = get_state_copy(pinfo.truth_file, pinfo.var, truth_index);
         plot3( ts(:,1), ts(:,2), ts(:,3), 'b'); hold on;
         legendstr = 'True State';
      end

      ens  = get_state_copy(pinfo.diagn_file, pinfo.var, ens_mean_index);
      plot3(ens(:,1), ens(:,2), ens(:,3), 'r');
      title(sprintf('%s Attractors for %s', pinfo.model, pinfo.diagn_file), ...    
                 'interpreter','none','fontweight','bold');

      if (exist('legendstr'))
         legend(legendstr,'Ensemble Mean',0);
      else
         legend(          'Ensemble Mean',0);
      end

      xlabel('state variable 1');
      ylabel('state variable 2');
      zlabel('state variable 3');
      legend boxoff

   case {'lorenz_96', 'lorenz_96_2scale', 'forced_lorenz_96', 'lorenz_04', ...
         'ikeda', 'simple_advection'}

      % Use one figure with subplots 
      figure(1); clf; iplot = 0;
      for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            subplot(length(pinfo.var_inds), 1, iplot);

            if (have_truth)
               truth    = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ivar);
               plot(pinfo.truth_times, truth,'b','LineWidth',2); hold on;
               legendstr = 'True State';
            end

            ens_mean    = get_var_series(pinfo.diagn_file, pinfo.var, ens_mean_index, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2));
            ens_members = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                         pinfo.diagn_time(1), pinfo.diagn_time(2));
            nmembers    = size(ens_members,2);

            plot(times,   ens_mean,'r','LineWidth',2); hold on;
            plot(times,ens_members,'g');

            if (exist('legendstr'))
               legend(legendstr, 'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
               plot(pinfo.truth_times,   truth,'b','LineWidth',2); % again, to put 'on top'
            else
               legend(           'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
            end
            plot(times,   ens_mean,'r','LineWidth',2.0); % again, to put 'on top'

            title(sprintf('%s %s varnum %d Ensemble Members of %s',...
                     pinfo.model, pinfo.var, ivar, pinfo.diagn_file), ...
                     'interpreter','none','fontweight','bold')
            xlabel(sprintf('model time (%d timesteps)',num_times));
            legend boxoff
      end

   case {'fms_bgrid','pe2lyr'}

      clf;

      timeunits = nc_attget(pinfo.fname, 'time',    'units');
      varunits  = nc_attget(pinfo.fname, pinfo.var, 'units');

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)

         if ( have_truth )
            truth    = GetCopy(pinfo.truth_file, truth_index,      pinfo );
            plot(pinfo.truth_times, truth,'b','LineWidth',2); hold on;
            legendstr = 'True State';
            s1 = sprintf('%s model ''%s'' Truth and %s Ensemble Members ', ...
                            pinfo.model, pinfo.var, pinfo.diagn_file);
         end

         ens_mean    = GetCopy(pinfo.diagn_file, ens_mean_index,   pinfo );
         ens_members = GetEns( pinfo.diagn_file,                   pinfo );
         nmembers    = size(ens_members,2);

         plot(times,   ens_mean,'r','LineWidth',2); hold on;
         plot(times,ens_members,'g');

         if (exist('legendstr'))
            legend(legendstr, 'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
            plot(pinfo.truth_times,   truth,'b','LineWidth',2); % again, to put 'on top'
         else
            legend(           'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
         end
         plot(times,   ens_mean,'r','LineWidth',2); % again, to put 'on top'

         s1 = sprintf('%s model ''%s'' %s Ensemble Members ', ...
                            pinfo.model, pinfo.var, pinfo.diagn_file);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                    pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1,s2},'interpreter','none','fontweight','bold');
         xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times));
         ylabel(varunits);
         legend boxoff

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

         subplot(2,1,2)

            if ( have_truth )
               truth     = GetCopy(pinfo.truth_file, truth_index,   pinfo );
               plot(pinfo.truth_times, truth,'b','LineWidth',2); hold on;
               legendstr = 'True State';
            end

            ens_mean    = GetCopy(pinfo.diagn_file, ens_mean_index, pinfo );
            ens_members = GetEns( pinfo.diagn_file,                 pinfo );
            nmembers    = size(ens_members,2);

            plot(times,   ens_mean,'r','LineWidth',2); hold on;
            plot(times,ens_members,'g'); 

            if (exist('legendstr'))
               legend(legendstr, 'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
               plot(pinfo.truth_times,   truth,'b','LineWidth',2); % again, to put 'on top'
            else
               legend(           'Ensemble Mean', sprintf('Ensemble Members (%d)',nmembers), 0);
            end

            s1 = sprintf('%s model ''%s'' %s Ensemble Members ', ...
                               pinfo.model, pinfo.var, pinfo.diagn_file);
            s2 = sprintf('level index %d lat %.2f lon %.2f', ...
                       pinfo.levelindex, pinfo.latitude, pinfo.longitude);
            title({s1,s2},'interpreter','none','fontweight','bold');
            xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times));
            ylabel(varunits);
            legend boxoff
            hold off;
      end

   otherwise

      error('model %s unknown.',pinfo.model)

end


%======================================================================
% Subfunctions
%======================================================================



function x = dim_length(fname,dimname)
y = nc_isvar(fname,dimname);
if (y < 1)
   error('%s has no %s dimension/coordinate variable',fname,varname)
end
bob = nc_getdiminfo(fname,dimname);
x   = bob.Length;



function var = GetCopy(fname, copyindex, pinfo)
% Gets a time-series of a single specified copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.levelindex = pinfo.levelindex;
myinfo.latindex   = pinfo.latindex;
myinfo.lonindex   = pinfo.lonindex;
[start, count]    = GetNCindices(myinfo,'diagn',pinfo.var);

var = nc_varget(fname, pinfo.var, start, count);






function var = GetEns(fname, pinfo)
% Gets a time-series of all copies of a prognostic variable 
% at a particular 3D location (level, lat, lon).
% Determining just the ensemble members (and not mean, spread ...)
% is the hard part.

% find which are actual ensemble members
metadata    = nc_varget(fname,'CopyMetaData');       % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's

if ( isempty(copyindices) )
   fprintf('%s has no valid ensemble members\n',fname)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',fname, num_copies)
   error('netcdf file has no ensemble members.')
end
ens_num     = length(copyindices);

myinfo.diagn_file = fname;
myinfo.levelindex = pinfo.levelindex;
myinfo.latindex   = pinfo.latindex;
myinfo.lonindex   = pinfo.lonindex;
[start, count]    = GetNCindices(myinfo,'diagn',pinfo.var);

bob = nc_varget(fname, pinfo.var, start, count); % 'bob' is only 2D 
var = bob(:,copyindices);



function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axis([0 360 -90 90])
   worldmap;
   axis image
   grid on



