function PlotJeffCorrel( pinfo )
%% Plots exploratory correlation plots. Don't use without talking to J. Anderson.
%
% Plots time series of correlations between two given variables.
% The correlation is done across ensemble members.
%
% PlotJeffCorrel is intended to be called by 'plot_jeff_correl'
%
% USAGE: PlotJeffCorrel(pinfo)
%
% pinfo is a model-dependent structure.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(pinfo.fname,'file') ~= 2), error('%s does not exist.',pinfo.fname), end

switch lower(pinfo.model)

   case {'fms_bgrid','pe2lyr','mitgcm_ocean','wrf','cam','sqg','pop'}
      clf;

      base_mem = get_hyperslab('fname',      pinfo.fname, ...
                               'varname',    pinfo.base_var, ...
                               'squeeze',    'T', ...
                               'permute',    'T', ...
                               'levelindex', pinfo.base_lvlind, ...
                               'latindex',   pinfo.base_latind, ...
                               'lonindex',   pinfo.base_lonind );

      comp_mem = get_hyperslab('fname',      pinfo.fname, ...
                               'varname',    pinfo.comp_var, ...
                               'squeeze',    'T', ...
                               'permute',    'T', ...
                               'levelindex', pinfo.comp_lvlind, ...
                               'latindex',   pinfo.comp_latind, ...
                               'lonindex',   pinfo.comp_lonind );

      correl = jeff_correl(base_mem, comp_mem);

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
         plot(pinfo.time, correl);

      s1 = sprintf('%s Correlation of ''%s'', lvl = %d, lat = %.2f, lon=%.2f', ...
          pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
          pinfo.base_lat, pinfo.base_lon);

      s2 = sprintf('with ''%s'', lvl = %d, lat = %.2f, lon= %.2f, %d ensemble members', ...
          pinfo.comp_var, pinfo.comp_lvl, pinfo.comp_lat, pinfo.comp_lon, ...
          pinfo.num_ens_members);

      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xdates(pinfo.time)
      ylabel('correlation')

   case {'mpas_atm'}
      error('not supported yet')

      clf;

      base_mem = get_hyperslab('fname',      pinfo.fname, ...
                               'varname',    pinfo.base_var, ...
                               'squeeze',    'T', ...
                               'permute',    'T', ...
                               'levelindex', pinfo.base_lvlind, ...
                               'cellindex',  pinfo.base_cellindex);

      comp_mem = get_hyperslab('fname',      pinfo.fname, ...
                               'varname',    pinfo.comp_var, ...
                               'squeeze',    'T', ...
                               'permute',    'T', ...
                               'levelindex', pinfo.comp_lvlind, ...
                               'cellindex',  pinfo.comp_cellindex);

      correl = jeff_correl(base_mem, comp_mem);

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
         plot(pinfo.time, correl);

      s1 = sprintf('%s Correlation of ''%s'', lvl = %d, lat = %.2f, lon=%.2f', ...
          pinfo.model, pinfo.base_var, pinfo.base_lvl, ...
          pinfo.base_lat, pinfo.base_lon);

      s2 = sprintf('with ''%s'', lvl = %d, lat = %.2f, lon= %.2f, %d ensemble members', ...
          pinfo.comp_var, pinfo.comp_lvl, pinfo.comp_lat, pinfo.comp_lon, ...
          pinfo.num_ens_members);

      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      datetick('x','yyyymmmdd HH:MM')
      ylabel('correlation')

   otherwise

      [num_vars,~] = nc_dim_info(pinfo.fname,'location');
      fprintf('PlotJeffCorrel: num_vars is %d\n',num_vars)

      % The Base Variable Index must be a valid state variable

      if ( pinfo.base_var_index > num_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, num_vars)
         error('you wanted variable # %d ', pinfo.base_var_index)
      end

      % The State Variable Index must also be a valid state variable

      if ( pinfo.base_var_index > num_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, num_vars)
         error('you wanted variable # %d ', pinfo.base_var_index)
      end

      % Get 'standard' ensemble series
      base_var  = get_hyperslab('fname',     pinfo.fname, ...
                               'varname',    pinfo.base_var, ...
                               'stateindex', pinfo.base_var_index, ...
                               'squeeze',    'T', ...
                               'permute',    'T');

      state_var = get_hyperslab('fname',     pinfo.fname, ...
                               'varname',    pinfo.state_var, ...
                               'stateindex', pinfo.state_var_index, ...
                               'squeeze',    'T', ...
                               'permute',    'T');

      % perform correlation
      correl = jeff_correl(base_var, state_var);

      clf; plot(pinfo.time,correl);

      s1 = sprintf('%s Correlation of variable %s %d, with variable %s %d', ...
               pinfo.model, pinfo.base_var, pinfo.base_var_index, ...
                      pinfo.state_var, pinfo.state_var_index);
      s2 = sprintf('%d ensemble members', pinfo.num_ens_members);
      title({s1,s2,pinfo.fname},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model "days" (%d timesteps)',pinfo.time_series_length))
      ylabel('correlation')

end


%======================================================================
% Subfunctions
%======================================================================


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
   datetick('x','mm.dd.HH','keeplimits'); % 'mm/dd'
   monstr = datestr(dates(1),31);
   xlabelstring = sprintf('month/day/HH - %s start',monstr);
else
   datetick('x',31,'keeplimits'); %'yyyy-mm-dd HH:MM:SS'
   monstr = datestr(dates(1),31);
   xlabelstring = sprintf('%s start',monstr);
end
xlabel(xlabelstring)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
