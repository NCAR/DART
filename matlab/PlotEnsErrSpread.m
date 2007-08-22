function PlotEnsErrSpread( pinfo )
% PlotEnsErrSpread     Creates summary plots of error and spread 
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

pinfo = CheckModelCompatibility(pinfo)

% Get some information from the truth_file 
ft = netcdf(pinfo.truth_file);
tmodel      = ft.model(:);
tvar_atts   = dim(ft{pinfo.var});      % cell array of dimensions for the var
%tnum_times  = length(tvar_atts{1});    % determine # of output times
tnum_copies = length(tvar_atts{2});    % # of ensemble members
tnum_vars   = length(tvar_atts{3});    % dimension of desired variable
close(ft);
% time might be a subset now of the entire array - compute it instead of
% assuming the entire size.
tnum_times = pinfo.truth_time(2) - pinfo.truth_time(1) + 1;

% Get some information from the diagn_file 
fd = netcdf(pinfo.diagn_file);
dmodel      = fd.model(:);
dvar_atts   = dim(fd{pinfo.var});      % cell array of dimensions for the var
%dnum_times  = length(dvar_atts{1});    % determine # of output times
dnum_copies = length(dvar_atts{2});    % # of ensemble members
dnum_vars   = length(dvar_atts{3});    % dimension of desired variable
close(fd);
% time might be a subset now of the entire array - compute it instead of
% assuming the entire size.
dnum_times = pinfo.diagn_time(2) - pinfo.diagn_time(1) + 1;

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Get some useful plotting arrays
times = getnc(pinfo.truth_file,'time', ...
              [ pinfo.truth_time(1) ], [ pinfo.truth_time(2) ]);

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

            disp(sprintf('%s model Variable %d',tmodel,ivar))

            subplot(3, 1, j);
               plot(times,err, 'b', ...
                    times,ens_spread(:, ivar), 'r');
               s1 = sprintf('%s model Var %d Ensemble Error Spread for %s', ...
                            tmodel, ivar, pinfo.diagn_file);
               title(s1,'interpreter','none','fontweight','bold')
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
               s1 = sprintf('%s model Var %d Ensemble Error Spread for %s', ...
                            tmodel, ivar, pinfo.diagn_file);
               title(s1,'interpreter','none','fontweight','bold')
               legend(string1,string2,0)
               legend boxoff
               xlabel(sprintf('model time (%d timesteps)',tnum_times))
               ylabel('distance')
      end

   case {'fms_bgrid','pe2lyr'}

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
         s1 = sprintf('%s model Var %s Ensemble Error Spread for %s', ...
                            tmodel, pinfo.var, pinfo.diagn_file);
         title(s1,'interpreter','none','fontweight','bold');

      title({ ...
        sprintf('Ensemble Mean Error, Ensemble Spread  %s ''%s'' for %s', ...
                tmodel, pinfo.var, pinfo.diagn_file), ...
        sprintf('level %d lat %.2f lon %.2f',pinfo.level, pinfo.latitude, ...
                 pinfo.longitude)}, 'interpreter','none','fontweight','bold');

         legend(string1,string2,0);
         legend boxoff
         xlabel(sprintf('model time (%d timesteps)',tnum_times));
         ylabel('distance');

   otherwise
      error(sprintf('model %s unknown.',tmodel))
end

%======================================================================
% Subfunctions
%======================================================================

function var = GetCopy(fname, copyindex, pinfo, tstart, tend)
% Gets a time-series of a single specified copy of a prognostic variable
% at a particular 3D location (level, lat, lon)
if strcmp(pinfo.var,'ps')
   corner = [tstart copyindex                  pinfo.latindex pinfo.lonindex];
   endpnt = [tend   copyindex                  pinfo.latindex pinfo.lonindex];
else
   corner = [tstart copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
   endpnt = [tend   copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
end
var = getnc(fname, pinfo.var, corner, endpnt);


function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axis([0 360 -90 90])
   worldmap;
   axis image
   grid on
