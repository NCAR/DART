function PlotCorrel( pinfo )
% PlotCorrel   space-time series of correlation between a variable at a given
% time and all variables at all times in an ensemble time sequence.
%
% PlotCorrel is intended to be called by 'plot_correl'.
%
% USAGE: PlotCorrel( pinfo )
%
% pinfo      A structure containing all necessary plotting information.
%            For the low-order models, the structure MUST contain:
%
% fname             name of netCDF file containing a DART ensemble
% base_var          name of netCDF variable
% base_var_index    index of state variable used as standard in correlation
% base_time         index of time series to use as the standard for correlation
%
% Example 1   (9var model with 1000 time steps)
%%------------------------------------------------------------------
% pinfo.fname          = 'Prior_Diag.nc';
% pinfo.base_var       = 'state';
% pinfo.base_var_index = 5;          % picked arbitrarily
% pinfo.base_time      = 238;        % ditto
% PlotCorrel(pinfo)                  % generates a plot

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

if (exist(pinfo.fname,'file') ~= 2), error('%s does not exist.',pinfo.fname), end

% Get some file-specific information.
model      = nc_attget(pinfo.fname,nc_global,'model');
timeunits  = nc_attget(pinfo.fname,'time','units');
varinfo    = nc_getvarinfo(pinfo.fname, pinfo.base_var);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}) )
      case 'time'
         num_times = varinfo.Size(i);
      case 'copy'
         num_copies = varinfo.Size(i);
      otherwise
         num_vars = varinfo.Size(i);
   end
end

switch(lower(model))

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection'}

      base_var_index = pinfo.base_var_index;
      base_time      = pinfo.base_time;
      
      % The Base Variable Index must be a valid state variable
      if ( base_var_index > num_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, num_vars)
         error('you wanted variable # %d ', base_var_index)
      end
      
      % The Time must be within range also.
      if ( base_time > num_times )
         fprintf('%s only has %d output times\n', pinfo.fname, num_times)
         error('you wanted time # %d ', base_time)
      end
      
      % Get 'standard' ensemble series 
      base = get_ens_series(pinfo.fname, pinfo.base_var, base_var_index);
      
      % It is efficient to preallocate correl storage ... 
      correl = zeros(num_vars,num_times);
      
      % Need to loop through all variables in the ensemble
      for i = 1:num_vars,
         state_var = get_ens_series(pinfo.fname, pinfo.base_var, i);
         correl(i, :) = ens_correl(base, base_time, state_var);
      end
     
      % Now for the plotting part ...
      disp('Please be patient ... this usually takes a bit ...')
      clf;
      
      contour(correl,-1:0.2:1);
      s1 = sprintf('%s Correlation of variable %s index %d, T = %d of %s', ...
               model, pinfo.base_var, base_var_index, base_time, pinfo.fname);
      s2 = sprintf('against all variables, all times, %d ensemble members', ...
               size(state_var,2)); 
      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel('time (timestep #)')
      ylabel('state variable (index)')
      set(gca,'YTick',1:num_vars)
      colorbar
      
      % highlight the reference state variable and time
      
      hold on;
      plot(base_time,base_var_index,'kh','MarkerSize',12,'MarkerFaceColor','k')

   case 'fms_bgrid'

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl   

      clf;

      times = nc_varget(pinfo.fname, 'time');
      switch lower(pinfo.comp_var)
         case {'ps','t'}
            lats     = nc_varget(pinfo.fname,'TmpJ'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'TmpI'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'TmpJ','units');
            lonunits = nc_attget(pinfo.fname,'TmpI','units');
         otherwise
            lats     = nc_varget(pinfo.fname,'VelJ'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'VelI'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'VelJ','units');
            lonunits = nc_attget(pinfo.fname,'VelI','units');
      end

      nxny = nx*ny;

      base_mem = Get1Ens( pinfo.fname, pinfo.base_var, pinfo.base_tmeind, ... 
                    pinfo.base_lvlind, pinfo.base_latind, pinfo.base_lonind );

      comp_ens = GetEnsLevel( pinfo.fname,       pinfo.comp_var, ...
                              pinfo.base_tmeind, pinfo.comp_lvlind);

      nmembers = size(comp_ens,1);

      corr = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end 

      correl = reshape(corr,[ny nx]);

      contour(lons,lats,correl,-1:0.2:1); hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %f of %s', ...
           model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, pinfo.base_time, pinfo.fname);

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, nmembers); 
      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      worldmap;
      axis image
      h = colorbar; 
      ax = get(h,'Position');
     %set(h,'Position',[ax(1) ax(2) ax(3)/2 ax(4)]);

   case 'mitgcm_ocean'

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl   

      clf;

      times      = nc_varget(pinfo.fname,'time');

      switch lower(pinfo.comp_var)
         case {'u'}
            lats     = nc_varget(pinfo.fname,'YC'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'XG'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'YC','units');
            lonunits = nc_attget(pinfo.fname,'XG','units');
         case {'v'}
            lats     = nc_varget(pinfo.fname,'YG'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'XC'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'YG','units');
            lonunits = nc_attget(pinfo.fname,'XC','units');
         otherwise
            lats     = nc_varget(pinfo.fname,'YC'); ny = length(lats);
            lons     = nc_varget(pinfo.fname,'XC'); nx = length(lons);
            latunits = nc_attget(pinfo.fname,'YC','units');
            lonunits = nc_attget(pinfo.fname,'XC','units');
      end

      nxny = nx*ny;

      base_mem = Get1Ens( pinfo.fname, pinfo.base_var, pinfo.base_tmeind, ... 
                    pinfo.base_lvlind, pinfo.base_latind, pinfo.base_lonind );
      comp_ens = GetEnsLevel( pinfo.fname,       pinfo.comp_var, ...
                              pinfo.base_tmeind, pinfo.comp_lvlind);
      nmembers = size(comp_ens,1);

      corr = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end 

      correl = reshape(corr,[ny nx]);

      contour(lons,lats,correl,-1:0.2:1); hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %f of %s', ...
           model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, pinfo.base_time, pinfo.fname);

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, nmembers); 
      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      worldmap;
      axis image
      h = colorbar; 
      ax = get(h,'Position');
     %set(h,'Position',[ax(1) ax(2) ax(3)/2 ax(4)]);

   case 'pe2lyr'

      % We are going to correlate one var/time/lvl/lat/lon  with
      % all other lats/lons for a var/time/lvl   

      clf;

      lats     = nc_varget(pinfo.fname,'lat'); ny = length(lats);
      lons     = nc_varget(pinfo.fname,'lon'); nx = length(lons);
      times    = nc_varget(pinfo.fname,'time');
      latunits = nc_attget(pinfo.fname,'lat','units');
      lonunits = nc_attget(pinfo.fname,'lon','units');

      nxny     = nx*ny;

      base_mem = Get1Ens( pinfo.fname, pinfo.base_var,    pinfo.base_tmeind, ... 
                    pinfo.base_lvlind, pinfo.base_latind, pinfo.base_lonind );
      comp_ens = GetEnsLevel( pinfo.fname,       pinfo.comp_var, ...
                              pinfo.base_tmeind, pinfo.comp_lvlind);
      nmembers = size(comp_ens,1);

      corr = zeros(nxny,1);

      for i = 1:nxny,
         x = corrcoef(base_mem, comp_ens(:, i));
         corr(i) = x(1, 2);
      end 

      correl = reshape(corr,[ny nx]);

      contour(lons,lats,correl,-1:0.2:1); hold on;
      plot(pinfo.base_lon, pinfo.base_lat, 'pk', ...
                 'MarkerSize',12,'MarkerFaceColor','k');
      s1 = sprintf('%s Correlation of ''%s'', level %d, (%.2f,%.2f) T = %f of %s', ...
           model, pinfo.base_var, pinfo.base_lvl, ...
             pinfo.base_lat, pinfo.base_lon, pinfo.base_time, pinfo.fname);

      s2 = sprintf('against ''%s'', entire level %d, same time, %d ensemble members', ...
               pinfo.comp_var, pinfo.comp_lvl, nmembers); 
      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel(sprintf('longitude (%s)',lonunits),'interpreter','none')
      ylabel(sprintf('latitude (%s)',latunits),'interpreter','none')
      worldmap;
      axis image
      h = colorbar; 
      % ax = get(h,'Position');
      % set(h,'Position',[ax(1) ax(2) ax(3)/2 ax(4)]);

   otherwise

      error('model %s not implemented yet', model)

end


%----------------------------------------------------------------------
% helper functions
%----------------------------------------------------------------------

function slice = Get1Ens(fname, var, tmeind, lvlind, latind, lonind)
% netcdf variable ordering is unimportant
% Get1Ens retrieves all the ensemble members for a particular 4D 
% location (time, level, lat, lon).
% The ensemble members do not include the mean, spread, etc. 

% find which are actual ensemble members
metadata    = nc_varget(fname,'CopyMetaData');           % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's

if ( isempty(copyindices) )
   fprintf('%s has no valid ensemble members\n',fname)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',fname, num_copies)
   error('netcdf file has no ensemble members.')
end

pinfo.diagn_file = fname;
pinfo.timeindex  = tmeind;
pinfo.levelindex = lvlind;
pinfo.latindex   = latind;
pinfo.lonindex   = lonind;
[start, count]   = GetNCindices(pinfo,'diagn',var);

% Get all the ensemble members and just subset the copies we want.

bob   = nc_varget(fname, var, start, count);
slice = bob(copyindices);



function slice = GetEnsLevel(fname, var, tmeind, lvlind)
% netcdf variable ordering is unimportant, nominally [ time copy level lat lon ]
% GetEnsLevel retrieves all the ensemble members for a particular time and level. 
%
% The ensemble members do not include the mean, spread, etc. 
% The level is returned as a [Nmem -by- NspatialLocations] matrix;
% each row is an observation and each column is a (spatial)variable. 
% should vectorize better.

% find which are actual ensemble members
metadata    = nc_varget(fname,'CopyMetaData');           % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's

if ( isempty(copyindices) )
   fprintf('%s has no valid ensemble members\n',fname)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',fname, num_copies)
   error('netcdf file has no ensemble members.')
end

pinfo.diagn_file = fname;
pinfo.timeindex  = tmeind;
pinfo.levelindex = lvlind;
[start, count]   = GetNCindices(pinfo,'diagn',var);

bob        = nc_varget(fname, var, start, count);
ted        = bob(copyindices,:,:);
[nm,ny,nx] = size(ted);
slice      = reshape(ted,[nm ny*nx]);



function PlotLocator(pinfo)
   plot(pinfo.base_lon, pinfo.base_lat,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axis([0 360 -90 90])
   worldmap
   axis image
   grid on
