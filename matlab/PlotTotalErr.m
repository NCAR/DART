function PlotTotalErr( pinfo )
% PlotTotalErr Plots summary plots of global error and spread
%
% The error is a measure of the distance between the 
% 'ensemble mean' and the 'true_state'. The distance of the 
% spread is also plotted.
%
% PlotTotalErr is intended to be called by 'plot_total_err'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotTotalErr( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
%
% Example 1   (Lorenz_63, Lorenz_96, 9var ...)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Posterior_Diag.nc';
% PlotTotalErr( pinfo )

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% $Source$
% $Revision$
% $Date$

CheckModelCompatibility(pinfo.truth_file, pinfo.diagn_file)

% Get some information from the truth_file 
ft = netcdf(pinfo.truth_file);
t.model      = ft.model(:);
timeunits    = ft{'time'}.units(:);
t.num_vars   = ncsize(ft('StateVariable')); % determine # of state variables
t.num_copies = ncsize(ft('copy')); % determine # of ensemble members
t.num_times  = ncsize(ft('time')); % determine # of output times
close(ft);

% Get some information from the diagn_file 
fd = netcdf(pinfo.diagn_file);
d.model      = fd.model(:);
d.num_vars   = ncsize(fd('StateVariable')); % determine # of state variables
d.num_copies = ncsize(fd('copy')); % determine # of ensemble members
d.num_times  = ncsize(fd('time')); % determine # of output times
close(fd);

% Get the netcdf variable indices for desired "copies"
% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Get some useful plotting arrays
times = getnc(pinfo.truth_file,'time');

switch lower(t.model)

   case 'MysteryModel'

      disp(sprintf('unknown model %s -- doing nothing',t.model))

   case 'fms_bgrid'

      BgridTotalError( pinfo )

   otherwise

      % Get the appropriate netcdf variables
      truth  = get_state_copy(pinfo.truth_file,     truth_index);
      ens    = get_state_copy(pinfo.diagn_file,  ens_mean_index);
      spread = get_state_copy(pinfo.diagn_file,ens_spread_index);

      % Also need to compute the spread; zero truth for this and
      % compute distance from 0
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/d.num_times;
      spreadTotal= sum(err_spread)/d.num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      clf;
      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      title(sprintf('%s Total Error over all %d variables for %s',...
                    t.model, d.num_vars, pinfo.diagn_file), ...
            'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',t.num_times))
      ylabel('Total Error')

end




function BgridTotalError( pinfo )
% netcdf has no state vector, it has prognostic variables.
% We are going to plot the total error (over a horizontal slice) 
% for each variable.
%

ft        = netcdf(pinfo.truth_file);
model     = ft.model(:); 
timeunits = ft{'time'}.units(:);
close(ft);

nvars = 4;

% Since the models are "compatible", get the info from either one.
tlons    = getnc(pinfo.truth_file,  'TmpI'); num_tlons  = length(tlons );
tlats    = getnc(pinfo.truth_file,  'TmpJ'); num_tlats  = length(tlats );
vlons    = getnc(pinfo.truth_file,  'VelI'); num_vlons  = length(vlons );
vlats    = getnc(pinfo.truth_file,  'VelJ'); num_vlats  = length(vlats );
levels   = getnc(pinfo.truth_file, 'level'); num_levels = length(levels);
times    = getnc(pinfo.truth_file,  'time'); num_times  = length(times );
ens_mems = getnc(pinfo.diagn_file,  'copy'); ens_size   = length(ens_mems);


% Try to coordinate "time" ... a poor attempt, needs refining
ens_times     = getnc(pinfo.diagn_file, 'time'); 
num_ens_times = length(ens_times);
if num_ens_times < num_times
   times     =     ens_times;
   num_times = num_ens_times;
end

% Initialize storage for error averaging
rms          = zeros(num_times, nvars, num_levels);
spread_final = zeros(num_times, nvars, num_levels);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state' ); 
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = SphereWeights(tlons, tlats);   % Temperature Grid
vwts = SphereWeights(vlons, vlats);   % Velocity    Grid

% Can we afford to get the whole thing at once ???
%----------------------------------------------------------------------
% surface pressure has only one level.
% GetPS   returns a   num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

disp('Processing surface pressure ...')

ivar   = 1;
ilevel = 1;

truth      = GetPS(pinfo.truth_file,      truth_index);
ens        = GetPS(pinfo.diagn_file,   ens_mean_index);
spread     = GetPS(pinfo.diagn_file, ens_spread_index);

err        = total_err(              truth,    ens, twts(:) );
err_spread = total_err(zeros(size(spread)), spread, twts(:) );

         rms(:,ivar,ilevel) = err;         % spatial mean
spread_final(:,ivar,ilevel) = err_spread;  % spatial mean

clear truth ens spread err err_spread  

%----------------------------------------------------------------------
% temperature ...  num_times x num_levels x num_lats x num_lons
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_levels,     % Loop through all levels

   disp(sprintf('Processing level %d of %d ...',ilevel,num_levels))

   %-------------------------------------------------------------------
   % temperature ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 2;

   truth  = GetLevel(pinfo.truth_file, ivar,      truth_index, ilevel);
   ens    = GetLevel(pinfo.diagn_file, ivar,   ens_mean_index, ilevel);
   spread = GetLevel(pinfo.diagn_file, ivar, ens_spread_index, ilevel);

   err        = total_err(              truth,    ens, twts(:));
   err_spread = total_err(zeros(size(spread)), spread, twts(:));

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;
  
   %-------------------------------------------------------------------
   % u ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 3;

   truth  = GetLevel(pinfo.truth_file, ivar,      truth_index, ilevel);
   ens    = GetLevel(pinfo.diagn_file, ivar,   ens_mean_index, ilevel);
   spread = GetLevel(pinfo.diagn_file, ivar, ens_spread_index, ilevel);

   err        = total_err(              truth,    ens, vwts(:));
   err_spread = total_err(zeros(size(spread)), spread, vwts(:));
   
            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;
  
   %-------------------------------------------------------------------
   % temperature ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 4;

   truth  = GetLevel(pinfo.truth_file, ivar,      truth_index, ilevel);
   ens    = GetLevel(pinfo.diagn_file, ivar,   ens_mean_index, ilevel);
   spread = GetLevel(pinfo.diagn_file, ivar, ens_spread_index, ilevel);

   err        = total_err(              truth,    ens, vwts(:));
   err_spread = total_err(zeros(size(spread)), spread, vwts(:));
   
            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;
  
end % End of level loop

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% Surface Pressure ... only one level
%----------------------------------------------------------------------
figure(1); clf;

      ft = netcdf(pinfo.truth_file);
      varunits = ft{'ps'}.units(:);
      close(ft);

      ivar = 1;
      h1 = plot(times,          rms(:, ivar, 1), '-'); hold on;
      h2 = plot(times, spread_final(:, ivar, 1), '--');
      s1 = sprintf('%s ''ps'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, ivar, 1)));
      s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))

%----------------------------------------------------------------------
% Temperature
%----------------------------------------------------------------------
figure(2); clf;
      ivar = 2;
      ft = netcdf(pinfo.truth_file);
      varunits = ft{'t'}.units(:);
      close(ft);

      h1 = plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
      h2 = plot(times, squeeze(spread_final(:, ivar, :)),'--');
      s1 = sprintf('%s ''temperature'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))

%----------------------------------------------------------------------
% U wind
%----------------------------------------------------------------------
figure(3); clf;
      ivar = 3;
      ft = netcdf(pinfo.truth_file);
      varunits = ft{'u'}.units(:);
      close(ft);

      h1 = plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
      h2 = plot(times, squeeze(spread_final(:, ivar, :)),'--');
      s1 = sprintf('%s ''U'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))

%----------------------------------------------------------------------
% V wind
%----------------------------------------------------------------------
figure(4); clf;
      ivar = 4;
      ft = netcdf(pinfo.truth_file);
      varunits = ft{'v'}.units(:);
      close(ft);

      h1 = plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
      h2 = plot(times, squeeze(spread_final(:, ivar, :)),'--');
      s1 = sprintf('%s ''V'' Ensemble Mean for %s', model,pinfo.diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))

%----------------------------------------------------------------------
% helper functions
%----------------------------------------------------------------------
function slice = GetPS(fname,copyindex);
corner     = [-1, copyindex, -1, -1];
endpnt     = [-1, copyindex, -1, -1];
ted        = getnc(fname,'ps',corner,endpnt);
[nt,ny,nx] = size(ted);
slice      = reshape(ted,[nt ny*nx]);


function slice = GetLevel(fname,ivar,copyindex,ilevel);
if ivar == 2 
   varstring = 't';
elseif ivar == 3 
   varstring = 'u';
elseif ivar == 4 
   varstring = 'v';
else
   error(sprintf(' variable id %d out of bounds',ivar))
end
corner     = [-1, copyindex, ilevel, -1, -1];
endpnt     = [-1, copyindex, ilevel, -1, -1];
ted        = getnc(fname,varstring,corner,endpnt);
[nt,ny,nx] = size(ted);
slice      = reshape(ted,[nt ny*nx]);


function wts = SphereWeights(lons, lats)
% SphereWeights creates weights based on area ...
%
% lons,lats must be 1D arrays (in degrees)

nlons = length(lons);
nlats = length(lats);

if ( prod(size(lons)) ~= nlons )
   disp('longitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end
if ( prod(size(lats)) ~= nlats )
   disp('latitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end

rads    = zeros(nlats,1);               % Ensure lats is a column vector,
rads(:) = pi*lats/180.0;                % and convert to radians.
wts     = cos( rads ) * ones(1,nlons);  % Results in a [nlat-x-nlon] matrix.
wts     = wts ./ sum(wts(:));           % Normalize to unity.
