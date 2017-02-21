function BgridTotalError( pinfo )
%% -------------------------------------------------------------------
% Plot the total area-weighted error for each variable.
%---------------------------------------------------------------------

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

for ivar=1:pinfo.num_state_vars,

   fprintf('Processing %s ...\n', pinfo.vars{ivar} )

   rmse     = zeros(pinfo.time_series_length,1);
   sprd     = zeros(pinfo.time_series_length,1);
   varunits = ncreadatt(pinfo.truth_file, pinfo.vars{ivar}, 'units');

   % determine what grid the variable lives on
   % determine the number of levels

   nlevels = 1;

   varinfo = ncinfo(pinfo.diagn_file, pinfo.vars{ivar});

   for idim = 1:length(varinfo.Dimensions),
      dimname   = varinfo.Dimensions(idim).Name;
      dimlength = varinfo.Dimensions(idim).Length;
      switch lower(dimname)
         case {'tmpj', 'velj'}
            latitudes   = ncread(pinfo.diagn_file, dimname);
         case {'tmpi', 'veli'}
            longitudes  = ncread(pinfo.diagn_file, dimname);
         case {'lev'}
            nlevels     = dimlength;
      end
   end

   % Calculate weights for area-averaging.
   weights = SphereWeights(latitudes, longitudes);

   for itime=1:pinfo.time_series_length,
       
       truth  = get_hyperslab('fname', pinfo.truth_file, ...
                    'varname', pinfo.vars{ivar}, ...
                    'permute', 'T', ...
                    'timeindex', pinfo.truth_time(1)+itime-1);

       ens    = get_hyperslab('fname', pinfo.diagn_file, ...
                    'varname', sprintf('%s_mean',pinfo.vars{ivar}), ...
                    'permute', 'T', ...
                    'timeindex', pinfo.diagn_time(1)+itime-1);

       varname = sprintf('%s_sd',pinfo.vars{ivar});
       [var_present, ~] = nc_var_exists(pinfo.diagn_file, varname);
       if (var_present > 0)
           spread = get_hyperslab('fname', pinfo.diagn_file, ...
                    'varname', varname, ...
                    'permute', 'T', ...
                    'timeindex', pinfo.diagn_time(1)+itime-1);
       else
           spread = zeros(size(ens));
       end

      %% Calculate the weighted mean squared error for each level.
      %  tensors come back [nlev,nlat,nlon] - or - [nlat,nlon]

      sqerr  = (truth - ens).^2;
      sqsprd =    spread    .^2;

      if (nlevels > 1) % take the mean over the first dimension
         sqerr  = squeeze(mean(sqerr ,1));
         sqsprd = squeeze(mean(sqsprd,1));
      end

      %% Create the (weighted) mean squared error

      ms_err    = sum(sqerr(:)  .* weights);
      ms_spread = sum(sqsprd(:) .* weights);

      %% Take the square root of the mean squared error
      rmse(itime) = sqrt(ms_err);
      sprd(itime) = sqrt(ms_spread);

   end % loop over time

   %-------------------------------------------------------------------
   % Each variable in its own figure window
   %-------------------------------------------------------------------
   figure(ivar); clf;
      plot(pinfo.time,rmse,'-', pinfo.time,sprd,'--')
      s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(rmse));
      s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(sprd));
      legend(s);
      legend boxoff;
      grid on;
      xdates(pinfo.time)
      ylabel(sprintf('global-area-weighted rmse (%s)',varunits))
      s1 = sprintf('%s %s Ensemble Mean', pinfo.model,pinfo.vars{ivar});
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end % loop around variables

clear truth ens spread err XY_spread




function weights = SphereWeights(lats,lons)
%% SphereWeights creates [nlat*nlon,1] matrix of weights based on latitude
%
% lats,lons must be 1D arrays (in degrees)

nlats = length(lats);
nlons = length(lons);

if ( numel(lats) ~= nlats )
   disp('latitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end
if ( numel(lons) ~= nlons )
   disp('longitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end

rads    = zeros(nlats,1);               % Ensure lats is a column vector,
rads(:) = pi*lats/180.0;                % and convert to radians.
wts     = cos( rads ) * ones(1,nlons);  % Results in a [nlat-x-nlon] matrix.
wts     = wts ./ sum(wts(:));           % Normalize to unity.
weights = wts(:);




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
