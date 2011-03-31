function PlotEnsErrSpread( pinfo )
%% PlotEnsErrSpread     Creates summary plots of error and spread 
%
% PlotEnsErrSpread is intended to be called by 'plot_ens_err_spread'.
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: EnsErrSpread( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models 
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
% var             name of netCDF variable of interest
% var_inds        indices of variables of interest 
%
% Example 0   (9var  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'state';
% pinfo.var_inds   = [ 1 2 3 4 5 6 7 8 9 ];
% PlotEnsErrSpread(pinfo)
%
% Example 1   (Lorenz_96  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'state';
% pinfo.var_inds   = [ 3 4 36 39 22 ];
% PlotEnsErrSpread(pinfo)
%
% Example 2   (Lorenz_96_2scale  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'X';
% pinfo.var_inds   = [ 3 18 27 ];
% PlotEnsErrSpread(pinfo)
%
% Example 3 (FMS BGrid model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotEnsErrSpread(pinfo)

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

pinfo = CheckModelCompatibility(pinfo);

% Get some information from the truth_file 
tmodel      = nc_attget( pinfo.truth_file,nc_global,'model');
tnum_times  = pinfo.truth_time(2) - pinfo.truth_time(1) + 1;

% Get some information from the diagn_file 
dmodel      = nc_attget( pinfo.diagn_file,nc_global,'model');
dnum_times  = pinfo.diagn_time(2) - pinfo.diagn_time(1) + 1;

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Get some useful plotting arrays
times = nc_varget(pinfo.truth_file,'time', pinfo.truth_time(1)-1 , tnum_times);

switch lower(tmodel)

   case '9var'

      % Get the appropriate copies
      truth      = get_state_copy(pinfo.truth_file, pinfo.var, truth_index, ...
                                  pinfo.truth_time(1), pinfo.truth_time(2)) ;
      ens_mean   = get_state_copy(pinfo.diagn_file, pinfo.var, ens_mean_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
      ens_spread = get_state_copy(pinfo.diagn_file, pinfo.var, ens_spread_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3

            ivar = (i - 1)*3 + j;

            err         = total_err(ens_mean(:,ivar) , truth(:,ivar));  
            errTotal    = sum(err)/dnum_times;
            spreadTotal = sum(ens_spread(:,ivar))/dnum_times;
            string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
            string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

            fprintf('%s model Variable %d\n',tmodel,ivar)

            subplot(3, 1, j);
               plot(times,err, 'b', ...
                    times,ens_spread(:, ivar), 'r');
               s1 = sprintf('%s model Var %d Ensemble Error Spread', tmodel, ivar);
               title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
               legend(string1,string2,0)
               legend boxoff
               xlabel(sprintf('model time (%d timesteps)',tnum_times))
               ylabel('distance')
         end
      end

   case {'lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection'} 
      % Get the appropriate copies
      truth      = get_state_copy(pinfo.truth_file, pinfo.var, truth_index, ...
                                  pinfo.truth_time(1), pinfo.truth_time(2)) ;
      ens_mean   = get_state_copy(pinfo.diagn_file, pinfo.var, ens_mean_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
      ens_spread = get_state_copy(pinfo.diagn_file, pinfo.var, ens_spread_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

      clf; iplot = 0;
      for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            err         = total_err(ens_mean(:,ivar) , truth(:,ivar));  
            errTotal    = sum(err)/dnum_times;
            spreadTotal = sum(ens_spread(:,ivar))/dnum_times;
            string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
            string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

            subplot(length(pinfo.var_inds), 1, iplot);
               plot(times,err, 'b', ...
                    times,ens_spread(:, ivar), 'r');
               s1 = sprintf('%s model Var %d Ensemble Error Spread', tmodel, ivar);
               title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
               legend(string1,string2,0)
               legend boxoff
               xlabel(sprintf('model time (%d timesteps)',tnum_times))
               ylabel('distance')
      end

   case {'fms_bgrid','pe2lyr','mitgcm_ocean'}

      clf;

      truth      = GetCopy(pinfo.truth_file, truth_index,      pinfo, ...
                           pinfo.truth_time(1), pinfo.truth_time(2)) ;
      ens_mean   = GetCopy(pinfo.diagn_file, ens_mean_index,   pinfo, ...
                           pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
      ens_spread = GetCopy(pinfo.diagn_file, ens_spread_index, pinfo, ...
                           pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

      subplot(2,1,1)
         PlotLocator(pinfo);

      subplot(2,1,2)
         err         = total_err(ens_mean, truth);  
         errTotal    = sum(err)/dnum_times;
         spreadTotal = sum(ens_spread)/dnum_times;
         string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
         string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

         plot(times,err, 'b', times,ens_spread, 'r');

         s1 = sprintf('Ensemble Mean Error, Ensemble Spread %s ''%s''',tmodel,pinfo.var);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                       pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1, s2, pinfo.diagn_file},'interpreter','none','fontweight','bold');

         legend(string1,string2,0);
         legend boxoff
         xlabel(sprintf('model time (%d timesteps)',tnum_times));
         ylabel('distance');

   otherwise
      error('model %s unknown.',tmodel)
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
   axis([0 360 -90 90])
   worldmap;
   axis image
   grid on
