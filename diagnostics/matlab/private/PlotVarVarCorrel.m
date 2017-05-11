function PlotVarVarCorrel( pinfo )
%% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.
% The correlation is done across ensemble members.
%
% PlotVarVarCorrel is intended to be called by 'plot_var_var_correl'
%
% USAGE: PlotVarVarCorrel( pinfo )
%
% pinfo is a structure with the following necessary components:
% fname
% base_var_index    index of one of the state variables to correlate
% base_time         index of the time of interest (timestep)
% state_var_index   index of the other state variable of interest.
%
% Example  (lorenz 63 model with 1000 timesteps)
%%--------------------------------------------------------
% pinfo.fname = 'analysis.nc';
% pinfo.base_var          = 'state';
% pinfo.base_var_index    = 2;
% pinfo.base_time         = 500;
% pinfo.state_var         = 'state';
% pinfo.state_var_index   = 1;
% PlotVarVarCorrel( pinfo )

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(pinfo.fname,'file') ~= 2), error('%s does not exist.',pinfo.fname), end

[num_times, ~, num_vars] = parse_varshape(pinfo.fname, pinfo.base_var);

switch lower(pinfo.model)

   case {'fms_bgrid','pe2lyr','wrf','cam','sqg','pop'}

      clf;

      base_mem = get_hyperslab('fname', pinfo.fname, ...
                     'varname', pinfo.base_var, ...
                     'squeeze', 'T',...
                     'permute', 'T', ...
                     'levelindex', pinfo.base_lvlind, ...
                     'latindex', pinfo.base_latind, ...
                     'lonindex', pinfo.base_lonind );

      comp_mem = get_hyperslab('fname', pinfo.fname, ...
                     'varname', pinfo.comp_var, ...
                     'squeeze', 'T',...
                     'permute', 'T', ...
                     'levelindex', pinfo.comp_lvlind, ...
                     'latindex', pinfo.comp_latind, ...
                     'lonindex', pinfo.comp_lonind );

      nmembers = size(comp_mem,2);

      correl = ens_correl(base_mem, pinfo.base_tmeind, comp_mem);

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
         plot(pinfo.time,correl);

      s1 = sprintf('%s Correlation of ''%s'', T = %d, lvl = %d, lat = %.2f, lon=%.2f', ...
          pinfo.model, pinfo.base_var, pinfo.base_time, pinfo.base_lvl, ...
          pinfo.base_lat, pinfo.base_lon);

      s2 = sprintf('with ''%s'', lvl = %d, lat = %.2f, lon= %.2f, %d ensemble members', ...
          pinfo.comp_var, pinfo.comp_lvl, pinfo.comp_lat, pinfo.comp_lon, nmembers);

      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xdates(pinfo.time)
      ylabel('correlation')

      % call out the time index in question, and put a corr==0 reference line.
      ax = axis;
      hold on;
      plot([pinfo.base_time pinfo.base_time],[ -1 1 ],'k:', ...
           [ax(1)         ax(2)],[  0 0 ],'k:')

   case {'mpas_atm'}

      clf;

      base_mem = get_hyperslab('fname', pinfo.fname, ...
                     'varname',pinfo.base_var, ...
                     'copy1',pinfo.ensemble_indices(1), ...
                     'copycount',pinfo.num_ens_members, ...
                     'levelindex',pinfo.base_lvlind, ...
                     'cellindex',pinfo.base_cellindex);

      comp_mem = get_hyperslab('fname', pinfo.fname, ...
                     'varname',pinfo.comp_var, ...
                     'copy1',pinfo.ensemble_indices(1), ...
                     'copycount',pinfo.num_ens_members, ...
                     'levelindex',pinfo.comp_lvlind, ...
                     'cellindex',pinfo.comp_cellindex);

      nmembers = size(comp_mem,2);

      correl = ens_correl(base_mem, pinfo.base_tmeind, comp_mem);

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
         plot(pinfo.time,correl);

      s1 = sprintf('%s Correlation of ''%s'', T = %d, lvl = %d, lat = %.2f, lon=%.2f', ...
          pinfo.model, pinfo.base_var, pinfo.base_time, pinfo.base_lvl, ...
          pinfo.base_lat, pinfo.base_lon);

      s2 = sprintf('with ''%s'', lvl = %d, lat = %.2f, lon= %.2f, %d ensemble members', ...
          pinfo.comp_var, pinfo.comp_lvl, pinfo.comp_lat, pinfo.comp_lon, nmembers);

      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xdates(pinfo.time)
      ylabel('correlation')

      % call out the time index in question, and put a corr==0 reference line.
      ax = axis;
      hold on;
      plot([pinfo.base_time pinfo.base_time],[ -1 1 ],'k:', ...
           [ax(1)                     ax(2)],[  0 0 ],'k:')

   otherwise

      fprintf('PlotVarVarCorrel: num_vars is %d\n',num_vars)

      % The Base Variable Index must be a valid state variable

      if ( pinfo.base_var_index > num_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, num_vars)
         error('you wanted variable # %d ', pinfo.base_var_index)
      end

      % The Time must be within range also.

      if ( pinfo.base_time > num_times )
         fprintf('%s only has %d output times\n', pinfo.fname, num_times)
         error('you wanted time # %d ', pinfo.base_time)
      end

      % The State Variable Index must also be a valid state variable

      if ( pinfo.base_var_index > num_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, num_vars)
         error('you wanted variable # %d ', pinfo.base_var_index)
      end

      % Get 'standard' ensemble series
      base_var  = get_hyperslab('fname', pinfo.fname,...
                     'squeeze', 'T',...
                     'permute', 'T', ...
                     'varname', pinfo.base_var,  ...
                     'stateindex', pinfo.base_var_index);

      state_var = get_hyperslab('fname', pinfo.fname,...
                     'squeeze', 'T',...
                     'permute', 'T', ...
                     'varname', pinfo.state_var, ...
                     'stateindex', pinfo.state_var_index);

      nmembers  = size(state_var,2);

      % perform a single correlation
      correl = ens_correl(base_var, pinfo.base_time, state_var);

      clf; plot(pinfo.time,correl);

      s1 = sprintf('%s Correlation of variable %s %d, timestep = %d, with variable %s %d', ...
               pinfo.model, pinfo.base_var, pinfo.base_var_index, pinfo.base_time, ...
                      pinfo.state_var, pinfo.state_var_index);
      s2 = sprintf('%d ensemble members -- %s', nmembers, pinfo.fname);
      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
      ylabel('correlation')

      % call out the time index in question, and put a corr==0 reference line.
      ax = axis;
      hold on;
      plot([pinfo.time(pinfo.base_time) pinfo.time(pinfo.base_time)],[ -1 1 ],'k:', ...
           [ax(1)         ax(2)],[  0 0 ],'k:')

      %axis(ax)
end


%======================================================================
% Subfunctions
%======================================================================

function [nT, nC, n3] = parse_varshape(fname,varname)

fileinfo  = ncinfo(fname);
nvars     = length(fileinfo.Variables);
isvar     = 0;

for i = 1:nvars
    candidate = fileinfo.Variables(i).Name;
    if (strcmp(varname,candidate))
        isvar = 1;
        break
    end
end

if (isvar < 1)
   error('%s has no variable named %s ',fname,varname)
end

nT = 0;
nC = 0;
n3 = 0;

varinfo = ncinfo(fname,varname);

% for i = 1:length(varinfo.Dimension)
for i = 1:3  % only want/need the first 3 dimensions.
   switch( lower(varinfo.Dimensions(i).Name))
      case 'time'
         nT = varinfo.Size(i);
      case 'member'
         nC = varinfo.Size(i);
      case 'location'
         n3 = varinfo.Size(i);
      otherwise
         n3 = varinfo.Size(i);
   end
end



function PlotLocator(pinfo)

   plot(pinfo.base_lon, pinfo.base_lat,'pb','MarkerSize',12,'MarkerFaceColor','b');
   hold on;
   plot(pinfo.comp_lon, pinfo.comp_lat,'pr','MarkerSize',12,'MarkerFaceColor','r');
   hold off;
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
if (length(get(gca,'XTick')) > 6)
   datetick('x','mm.dd.HH','keeplimits','keepticks'); % 'mm/dd'
   monstr = datestr(dates(1),31);
   xlabelstring = sprintf('month/day/HH - %s start',monstr);
else
   datetick('x',31,'keeplimits','keepticks'); %'yyyy-mm-dd HH:MM:SS'
   monstr = datestr(dates(1),31);
   xlabelstring = sprintf('%s start',monstr);
end
xlabel(xlabelstring)


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
