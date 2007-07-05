function PlotBins(pinfo)
% PlotBins Plots rank histograms of ensemble mean
%
% PlotBins is intended to be called by 'plot_bins'.
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotBins(pinfo)
% 
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copy tagged 'ensemble mean'
% state_var_inds  indices of state variables of interest
%
% Example 1 (Lorenz_96  model)
%%-------------------------------------------------------- 
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.state_var_inds = [3 4 36 39 22];
% PlotBins( pinfo );
%
% Example 2 (FMS BGrid model)
%%-------------------------------------------------------- 
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotBins( pinfo );

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

% Get the state for the truth
truth_index = get_copy_index(pinfo.truth_file,'true state');
true_model  =         GetAtt(pinfo.truth_file,'model');

switch lower(true_model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3
            ivar = (i - 1)*3 + j;
            truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                             ivar, pinfo.truth_time(1), pinfo.truth_time(2));
            ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
            bins  = rank_hist(ens, truth);
            subplot(3, 1, j);
            bar(bins);
            title(sprintf('%s Variable %d for %s', ...
                  true_model,ivar,pinfo.diagn_file), ...
                  'interpreter','none','fontweight','bold')
         end
      end

   case {'lorenz_63','lorenz_84','lorenz_96','lorenz_04','forced_lorenz_96','ikeda'}

      clf; iplot = 0;
      for ivar = pinfo.state_var_inds,
         iplot = iplot + 1;
         truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                               ivar, pinfo.truth_time(1), pinfo.truth_time(2));
         ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
         bins  = rank_hist(ens, truth);
         subplot(length(pinfo.state_var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %d for %s', ...
               true_model,ivar,pinfo.diagn_file), ...
               'interpreter','none','fontweight','bold')
      end

   case {'lorenz_96_2scale','simple_advection'}

      clf; iplot = 0;
      for ivar = pinfo.state_var_inds,
         iplot = iplot + 1;

         truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                               ivar, pinfo.truth_time(1), pinfo.truth_time(2));
         ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
         bins  = rank_hist(ens, truth);
         subplot(length(pinfo.state_var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %s %d for %s', ...
               true_model,pinfo.var, ivar,pinfo.diagn_file), ...
               'interpreter','none','fontweight','bold')
         ax = axis;
         ax(1) = 0.5;
         ax(2) = length(bins)+0.5;
         axis(ax)
      end

   case 'fms_bgrid'

      % need to know which prognostic variable  1=ps 2=t 3=u 4=v
      % which level 1<=ps<=1     1<= t,u,v <= nlevel
      % which location

      clf;

      truth = GetCopy(pinfo.truth_file, truth_index, pinfo, ...
                      pinfo.truth_time(1), pinfo.truth_time(2));
      ens   = GetEns( pinfo.diagn_file, pinfo, ...
                      pinfo.diagn_time(1), pinfo.diagn_time(2));

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
      bins  = rank_hist(ens, truth);
      bar(bins);
      title({ ...
        sprintf('%s ''%s'' for %s ', true_model, pinfo.var, pinfo.diagn_file), ...
        sprintf('level %d lat %.2f lon %.2f',pinfo.level, pinfo.latitude, ...
                 pinfo.longitude)}, 'interpreter','none','fontweight','bold')

   otherwise

      error(sprintf('model %s unknown',model))

end


%======================================================================
% Subfunctions
%======================================================================

function modelstring = GetAtt(fname,attname)
% Get a global attribute from a netCDF file.
f = netcdf(fname,'nowrite');   % open with low-level netcdf operators.
modelstring = f.model(:);      % grab a global attribute
close(f)
if isempty(modelstring)
   error(sprintf('NO model in netCDF file %s',fname))
else
   disp(sprintf('Selected model is %s',modelstring))
end



function var = GetCopy(fname, copyindex, pinfo, tstart, tend)
% Gets a time-series of a single specified copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)

switch(lower(pinfo.var))
   case {'ps'}
      corner = [tstart copyindex                 pinfo.latindex pinfo.lonindex];
      endpnt = [tend   copyindex                 pinfo.latindex pinfo.lonindex];
   otherwise
      corner = [tstart copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
      endpnt = [tend   copyindex pinfo.levelindex pinfo.latindex pinfo.lonindex];
end
var = getnc(fname, pinfo.var, corner, endpnt);


function var = get_1Dvar_type_series(fname, copyindex, vrbl, vrbl_ind)
% Gets a time series of a single specified copy of a prognostic variable 
% The (spatially-) 1D vars are   (time,copy,location) 

% this seems unused, but if it is called, tstart and tend should be passed
% in and used instead of the -1 in the corner and endpt below.
corner = [-1 copyindex vrbl_ind ];
endpnt = [-1 copyindex vrbl_ind ];
var = getnc(fname, pinfo.var, corner, endpnt);




function var = GetEns(fname, pinfo, tstart, tend)
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
   corner = [tstart -1                  pinfo.latindex pinfo.lonindex];
   endpnt = [tend   -1                  pinfo.latindex pinfo.lonindex];
else
   corner = [tstart -1 pinfo.levelindex pinfo.latindex pinfo.lonindex];
   endpnt = [tend   -1 pinfo.levelindex pinfo.latindex pinfo.lonindex];
end
bob = getnc(fname, pinfo.var, corner, endpnt); % 'bob' is only 2D 
var = bob(:,copyindices);


function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pg','MarkerSize',12,'MarkerFaceColor','g');
   axis([0 360 -90 90])
   worldmap
   
