function PlotCEnsErrSpread( pinfo )
%% PlotCEnsErrSpread

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% this sets start/stop time indices for both truth and diagn file now
pinfo = CheckModelCompatibility(pinfo.truth_file, pinfo.diagn_file)

model     = ncreadatt(pinfo.truth_file, '/', 'model');
timeunits = ncreadatt(pinfo.truth_file,'time','units');

nvars = 4;

% Since the models are "compatible", get the info from either one.
levels    = ncread(pinfo.truth_file, 'level'); num_levels = length(levels);
ens_mems  = ncread(pinfo.diagn_file,  'copy'); ens_size   = length(ens_mems);
num_times = pinfo.time_series_length;

% Initialize storage for error averaging
rms      = zeros(num_times, nvars, num_levels);
sd_final = zeros(num_times, nvars, num_levels);

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_member_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_member_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_member_index(pinfo.diagn_file, 'ensemble spread');


% Can we afford to get the whole thing at once ???
%----------------------------------------------------------------------
% surface pressure has only one level.
% GetPS   returns a   num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

disp('Processing surface pressure ...')

ivar   = 1;
ilevel = 1;

field  = GetPS(pinfo.truth_file,      truth_index, ...
               pinfo.truth_time(1), pinfo.truth_time(2)) ;
ens    = GetPS(pinfo.diagn_file,   ens_mean_index, ...
               pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
sd     = GetPS(pinfo.diagn_file, ens_spread_index, ...
               pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

% Compute statistics of the UNWEIGHTED error field
% should weight by the area ...
ens_err = ens - field;

     rms(:,ivar,ilevel) = mean(abs(ens_err),2);   % spatial mean
sd_final(:,ivar,ilevel) = mean(sd,2);             % spatial mean

clear field ens sd ens_err

%----------------------------------------------------------------------
% temperature ...  num_times x num_levels x num_lats x num_lons
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_levels,     % Loop through all levels

   fprintf('Processing level %d of %d ...\n',ilevel,num_levels)

   %-------------------------------------------------------------------
   % temperature ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 2;

   field      = GetLevel(pinfo.truth_file, ivar,      truth_index, ilevel, ...
                         pinfo.truth_time(1), pinfo.truth_time(2)) ;
   ens        = GetLevel(pinfo.diagn_file, ivar,   ens_mean_index, ilevel, ...
                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
   sd         = GetLevel(pinfo.diagn_file, ivar, ens_spread_index, ilevel, ...
                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

   % Compute statistics of the UNWEIGHTED error field
   % should weight by the area ...
   ens_err = ens - field;

        rms(:, ivar, ilevel) = mean(abs(ens_err),2);
   sd_final(:, ivar, ilevel) = mean(sd,2);

   %-------------------------------------------------------------------
   % u ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 3;

   field      = GetLevel(pinfo.truth_file, ivar,      truth_index, ilevel, ...
                         pinfo.truth_time(1), pinfo.truth_time(2)) ;
   ens        = GetLevel(pinfo.diagn_file, ivar,   ens_mean_index, ilevel, ...
                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
   sd         = GetLevel(pinfo.diagn_file, ivar, ens_spread_index, ilevel, ...
                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

   % Compute statistics of the UNWEIGHTED error field
   % should weight by the area ...
   ens_err = ens - field;

        rms(:, ivar, ilevel) = mean(abs(ens_err),2);
   sd_final(:, ivar, ilevel) = mean(sd,2);

   %-------------------------------------------------------------------
   % temperature ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 4;

   field      = GetLevel(pinfo.truth_file, ivar,      truth_index, ilevel, ...
                         pinfo.truth_time(1), pinfo.truth_time(2)) ;
   ens        = GetLevel(pinfo.diagn_file, ivar,   ens_mean_index, ilevel, ...
                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
   sd         = GetLevel(pinfo.diagn_file, ivar, ens_spread_index, ilevel, ...
                         pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

   % Compute statistics of the UNWEIGHTED error field
   % should weight by the area ...
   ens_err = ens - field;

        rms(:, ivar, ilevel) = mean(abs(ens_err),2);
   sd_final(:, ivar, ilevel) = mean(sd,2);

end % End of level loop

clear field ens sd ens_err

%----------------------------------------------------------------------
% Surface Pressure ... only one level
%----------------------------------------------------------------------
figure(1); clf;
      ivar = 1;
      plot(pinfo.time, rms(:, ivar, 1), '-', pinfo.time, sd_final(:, ivar, 1), '--');

      s1 = sprintf('%s model ''ps'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      s1 = sprintf('time-mean ensemble mean error  = %f', mean(     rms(:, ivar, 1)));
      s2 = sprintf('time-mean ensemble spread = %f', mean(sd_final(:, ivar, 1)));
      h = legend(s1,s2); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel('distance')

%----------------------------------------------------------------------
% Temperature
%----------------------------------------------------------------------
figure(2); clf;
      ivar = 2;
      plot_temp = reshape(rms(:, ivar, :), [num_times num_levels]);

      h1 = plot(pinfo.time, squeeze(     rms(:, ivar, :)),'-'); hold on;
      h2 = plot(pinfo.time, squeeze(sd_final(:, ivar, :)),'--');
      s1 = sprintf('%s model ''temperature'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(     rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f',i,mean(sd_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel('distance')

%----------------------------------------------------------------------
% U wind
%----------------------------------------------------------------------
figure(3); clf;
      ivar = 3;
      plot_temp = reshape(rms(:, ivar, :), [num_times num_levels]);

      h1 = plot(pinfo.time, squeeze(     rms(:, ivar, :)),'-'); hold on;
      h2 = plot(pinfo.time, squeeze(sd_final(:, ivar, :)),'--');
      s1 = sprintf('%s model ''U'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(     rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f',i,mean(sd_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel('distance')

%----------------------------------------------------------------------
% V wind
%----------------------------------------------------------------------
figure(4); clf;
      ivar = 4;
      plot_temp = reshape(rms(:, ivar, :), [num_times num_levels]);

      h1 = plot(pinfo.time, squeeze(     rms(:, ivar, :)),'-'); hold on;
      h2 = plot(pinfo.time, squeeze(sd_final(:, ivar, :)),'--');
      s1 = sprintf('%s model ''V'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(     rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f',i,mean(sd_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel('distance')

%----------------------------------------------------------------------
% helper functions
%----------------------------------------------------------------------
function slice = GetPS(fname,copyindex,tstartind,tendind);

varinfo = nc_getvarinfo(fname,'ps');

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tendind - tstartind + 1;
      case{'copy'}
         start(i) = copyindex - 1;
         count(i) = 1;
      otherwise
         start(i) =  0;
         count(i) = -1;
   end
end

error('should use get_hyperslab() instead')
ted        = nc_varget(fname,'ps',start,count);
[nt,ny,nx] = size(ted);
slice      = reshape(ted,[nt ny*nx]);





function slice = GetLevel(fname,ivar,copyindex,ilevel,tstartind,tendind);
if ivar == 2
   varstring = 't';
elseif ivar == 3
   varstring = 'u';
elseif ivar == 4
   varstring = 'v';
else
   error(' variable id %d out of bounds',ivar)
end

varinfo = nc_getvarinfo(fname,varstring);
ndims   = length(varinfo.Dimension);
start   = zeros(1,ndims);
count   = zeros(1,ndims);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tendind - tstartind + 1;
      case{'copy'}
         start(i) = copyindex - 1;
         count(i) = 1;
      case{'lev'}
         start(i) = ilevel - 1;
         count(i) = 1;
      otherwise
         start(i) =  0;
         count(i) = -1;
   end
end
error('should use get_hyperslab() instead')
ted        = nc_varget(fname,varstring,start,count);
[nt,ny,nx] = size(ted);
slice      = reshape(ted,[nt ny*nx]);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
