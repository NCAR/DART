function PlotCEnsErrSpread(truth_file, diagn_file)
%
%
%

CheckModelCompatibility(truth_file,diagn_file)

ft     = netcdf(truth_file);
model  = ft.model(:); 
close(ft);

nvars = 4;

% Since the models are "compatible", get the info from either one.

tlons    = getnc(truth_file,  'TmpI'); num_tlons  = length(tlons );
tlats    = getnc(truth_file,  'TmpJ'); num_tlats  = length(tlats );
vlons    = getnc(truth_file,  'VelI'); num_vlons  = length(vlons );
vlats    = getnc(truth_file,  'VelJ'); num_vlats  = length(vlats );
levels   = getnc(truth_file, 'level'); num_levels = length(levels);
times    = getnc(truth_file,  'time'); num_times  = length(times );
ens_mems = getnc(diagn_file,  'copy'); ens_size   = length(ens_mems);

% Try to coordinate "time" ... a poor attempt, needs refining
ens_times     = getnc(diagn_file, 'time'); 
num_ens_times = length(ens_times);
if num_ens_times < num_times
   times     =     ens_times;
   num_times = num_ens_times;
end

% Initialize storage for error averaging
rms      = zeros(num_times, nvars, num_levels);
sd_final = zeros(num_times, nvars, num_levels);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(truth_file, 'true state' ); 
ens_mean_index   = get_copy_index(diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(diagn_file, 'ensemble spread');


% Can we afford to get the whole thing at once ???
%----------------------------------------------------------------------
% surface pressure has only one level.
% GetPS   returns a   num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

disp('Processing surface pressure ...')

ivar   = 1;
ilevel = 1;

field  = GetPS(truth_file,      truth_index);
ens    = GetPS(diagn_file,   ens_mean_index);
sd     = GetPS(diagn_file, ens_spread_index);

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

   disp(sprintf('Processing level %d of %d ...',ilevel,num_levels))

   %-------------------------------------------------------------------
   % temperature ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 2;

   field      = GetLevel(truth_file, ivar,      truth_index, ilevel);
   ens        = GetLevel(diagn_file, ivar,   ens_mean_index, ilevel);
   sd         = GetLevel(diagn_file, ivar, ens_spread_index, ilevel);

   % Compute statistics of the UNWEIGHTED error field
   % should weight by the area ...
   ens_err = ens - field;

        rms(:, ivar, ilevel) = mean(abs(ens_err),2);
   sd_final(:, ivar, ilevel) = mean(sd,2);
  
   %-------------------------------------------------------------------
   % u ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 3;

   field      = GetLevel(truth_file, ivar,      truth_index, ilevel);
   ens        = GetLevel(diagn_file, ivar,   ens_mean_index, ilevel);
   sd         = GetLevel(diagn_file, ivar, ens_spread_index, ilevel);

   % Compute statistics of the UNWEIGHTED error field
   % should weight by the area ...
   ens_err = ens - field;

        rms(:, ivar, ilevel) = mean(abs(ens_err),2);
   sd_final(:, ivar, ilevel) = mean(sd,2);
  
   %-------------------------------------------------------------------
   % temperature ...  num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------
   ivar   = 4;

   field      = GetLevel(truth_file, ivar,      truth_index, ilevel);
   ens        = GetLevel(diagn_file, ivar,   ens_mean_index, ilevel);
   sd         = GetLevel(diagn_file, ivar, ens_spread_index, ilevel);

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
      plot(times, rms(:, ivar, 1), '-', times, sd_final(:, ivar, 1), '--');

      s1 = sprintf('%s model ''ps'' Ensemble Mean for %s', model,diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      s1 = sprintf('time-mean ensemble mean error  = %f', mean(     rms(:, ivar, 1)));
      s2 = sprintf('time-mean ensemble spread = %f', mean(sd_final(:, ivar, 1)));
      h = legend(s1,s2); legend(h,'boxoff')
      grid on;

%----------------------------------------------------------------------
% Temperature
%----------------------------------------------------------------------
figure(2); clf;
      ivar = 2;
      plot_temp = reshape(rms(:, ivar, :), [num_times num_levels]);

      h1 = plot(times, squeeze(     rms(:, ivar, :)),'-'); hold on;
      h2 = plot(times, squeeze(sd_final(:, ivar, :)),'--');
      s1 = sprintf('%s model ''temperature'' Ensemble Mean for %s', model,diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(     rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f',i,mean(sd_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;

%----------------------------------------------------------------------
% U wind
%----------------------------------------------------------------------
figure(3); clf;
      ivar = 3;
      plot_temp = reshape(rms(:, ivar, :), [num_times num_levels]);

      h1 = plot(times, squeeze(     rms(:, ivar, :)),'-'); hold on;
      h2 = plot(times, squeeze(sd_final(:, ivar, :)),'--');
      s1 = sprintf('%s model ''U'' Ensemble Mean for %s', model,diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(     rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f',i,mean(sd_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;

%----------------------------------------------------------------------
% V wind
%----------------------------------------------------------------------
figure(4); clf;
      ivar = 4;
      plot_temp = reshape(rms(:, ivar, :), [num_times num_levels]);

      h1 = plot(times, squeeze(     rms(:, ivar, :)),'-'); hold on;
      h2 = plot(times, squeeze(sd_final(:, ivar, :)),'--');
      s1 = sprintf('%s model ''V'' Ensemble Mean for %s', model,diagn_file);
      title(s1,'interpreter','none','fontweight','bold')

      for i = 1:num_levels,
         s{i           } = sprintf('level %d error  %.3f', i,mean(     rms(:, ivar, i)));
         s{i+num_levels} = sprintf('level %d spread %.3f',i,mean(sd_final(:, ivar, i)));
      end
      h = legend([h1 h2],s); legend(h,'boxoff')
      grid on;

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


