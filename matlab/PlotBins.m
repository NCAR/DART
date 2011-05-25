function PlotBins(pinfo)
%% PlotBins Plots rank histograms of ensemble mean
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
% var_inds        indices of state variables of interest
%
% Example 1 (Lorenz_96  model)
%%-------------------------------------------------------- 
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var_inds   = [3 4 36 39 22];
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

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

pinfo = CheckModelCompatibility(pinfo)

% Get the state for the truth
truth_index = get_copy_index(pinfo.truth_file,'true state');
true_model  = nc_attget(pinfo.truth_file, nc_global, 'model');

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
            xlabel('rank')
            ylabel('occurrence')
            ax = axis;
            ax(1) = 0.5;
            ax(2) = length(bins)+0.5;
            axis(ax)
            axis tight
         end
      end

   case {'lorenz_63','lorenz_84','lorenz_96','lorenz_04','forced_lorenz_96','ikeda'}

      clf; iplot = 0;
      for ivar = pinfo.var_inds,
         iplot = iplot + 1;
         truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                               ivar, pinfo.truth_time(1), pinfo.truth_time(2));
         ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
         bins  = rank_hist(ens, truth);
         subplot(length(pinfo.var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %d for %s', ...
               true_model,ivar,pinfo.diagn_file), ...
               'interpreter','none','fontweight','bold')
         xlabel('rank')
         ylabel('occurrence')
         ax = axis;
         ax(1) = 0.5;
         ax(2) = length(bins)+0.5;
         axis(ax)
         axis tight
      end

   case {'lorenz_96_2scale','simple_advection'}

      clf; iplot = 0;
      for ivar = pinfo.var_inds,
         iplot = iplot + 1;

         truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                               ivar, pinfo.truth_time(1), pinfo.truth_time(2));
         ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
         bins  = rank_hist(ens, truth);
         subplot(length(pinfo.var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %s %d for %s', ...
               true_model,pinfo.var, ivar,pinfo.diagn_file), ...
               'interpreter','none','fontweight','bold')
         xlabel('rank')
         ylabel('occurrence')
         ax = axis;
         ax(1) = 0.5;
         ax(2) = length(bins)+0.5;
         axis(ax)
         axis tight
      end

   case {'fms_bgrid','pe2lyr','mitgcm_ocean'}

      % It is intended that all 3D models have all the required information
      % set in the corresponding Get<model>Info.m script.

      clf;

      truth = GetCopy(pinfo);
      ens   = GetEns( pinfo);

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
      bins  = rank_hist(ens, truth);
      bar(bins);
      title({ ...
        sprintf('%s ''%s'' for %s ', true_model, pinfo.var, pinfo.diagn_file), ...
        sprintf('level %d lat %.2f lon %.2f',pinfo.level, pinfo.latitude, ...
                 pinfo.longitude)}, 'interpreter','none','fontweight','bold')
      xlabel('rank')
      ylabel('occurrence')
      ax = axis;
      ax(1) = 0.5;
      ax(2) = length(bins)+0.5;
      axis(ax)
      axis tight

   otherwise

      error('model %s unknown',true_model)

end


%======================================================================
% Subfunctions
%======================================================================


function var = GetCopy(pinfo)
% Gets a time-series of a single specified 'true' copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)

[start, count] = GetNCindices(pinfo,'truth',pinfo.var);
varinfo        = nc_getvarinfo(pinfo.truth_file, pinfo.var);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = pinfo.truth_time(1) - 1;
         count(i) = pinfo.truth_time(2) - pinfo.truth_time(1) + 1;
         break
      otherwise
   end
end

var = nc_varget(pinfo.truth_file, pinfo.var, start, count);



function var = GetEns(pinfo)
% Gets a time-series of all copies of a prognostic variable 
% at a particular 3D location (level, lat, lon).
% Determining just the ensemble members (and not mean, spread ...)
% is the hard part.

% find which are actual ensemble members
metadata    = nc_varget(pinfo.diagn_file,'CopyMetaData');       % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's

if ( isempty(copyindices) )
   fprintf('%s has no valid ensemble members\n',pinfo.diagn_file)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',pinfo.diagn_file, size(metadata,1))
   error('netcdf file has no ensemble members.')
end
ens_num     = length(copyindices);

% Get all ensemble members, just return desired ones.
% This makes fewer assumptions about variable shape.
[start, count] = GetNCindices(pinfo,'diagn',pinfo.var);
varinfo        = nc_getvarinfo(pinfo.diagn_file, pinfo.var);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = pinfo.diagn_time(1) - 1;
         count(i) = pinfo.diagn_time(2) - pinfo.diagn_time(1) + 1;
         break
      otherwise
   end
end

bob = nc_varget(pinfo.diagn_file, pinfo.var, start, count); % 'bob' is only 2D time-X-copy
var = bob(:,copyindices);


% could/should check input for valid range, etc.


function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axis([0 360 -90 90]);
   worldmap;
   
