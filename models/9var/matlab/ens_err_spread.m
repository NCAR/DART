function ens_err_spread(dir_name,tname,ename)
% ens_err_spread
%
% can be called with no arguments (in which case you are queried)
% or can be called with exactly three arguments:
%
% ens_error_spread(DirName,TruthName,EnsName)
%
% DirName    directory containing data files (may be '.')
% TruthName  file containing the "True State" of the system.
% EnsName    file containing the "Ensemble" estimate of the system.
%
% Example:
%
% dirname  = '../work';
% truename = 'True_State.nc';
% ensname  = 'Prior_Diag.nc';
% ens_err_spread(dirname,truename,ensname);


if nargin < 0

   disp('\n Remember to enclose your responses in single quotes.\n')

   dir_name = input('Input directory; . for current directory  ');
   tname    = input('Input file name for True state  ');
   ename    = input('Input file name for ensemble  ');

end

%----------------------------------------------------------------------
% Load the "true state" file
%----------------------------------------------------------------------
fname = fullfile(dir_name,tname);

loc            = getnc(fname, 'loc1d');
num_loc        = size(loc, 1);
true_times     = getnc(fname, 'time');
num_true_times = size(true_times, 1);

state_vec      = getnc(fname, 'state');
state_size     = size(state_vec, 2);

%----------------------------------------------------------------------
% Load the ensemble "state" vector.  state(time, copy, StateVariable) 
% Assumes two copies are ensemble mean followed by ensemble spread
% Should be automated and checked at some point
%----------------------------------------------------------------------
ens_fname = fullfile(dir_name,ename);

ens_full_vec = getnc(ens_fname, 'state');    % get entire "state" vector
ens_size     = size(ens_full_vec, 2);        % Ensemble size 

if (ens_size < 2) 
   error(sprintf('%s does not contain an ensemble mean and an ensemble spread.',ens_fname))
end

ens_times     = getnc(ens_fname, 'time');       
num_ens_times = size(ens_times, 1);

%----------------------------------------------------------------------
% assumes both start at zero ... 

num_times = min(num_true_times, num_ens_times);   % can be different.

% Initialize storage for error averaging
rms(1:num_times) = 0.0;
sd_final(1:num_times) = 0.0;


% Loop through all the time levels
for time_ind = 1 : num_times

   % Extract state and ensemble for just this time
   single_state =    state_vec(time_ind,    :);
   ens_mean     = ens_full_vec(time_ind, 1, :);
   ens_sd       = ens_full_vec(time_ind, 2, :);

   % Compute the difference field
   ens_err = (ens_mean(:) - single_state(:));   % ens_err is a simple array

   % Compute statistics of the entire error field
   max_err = max(ens_err);
   min_err = min(ens_err);
   rms_err = mean(abs(ens_err));
   sd_mean = mean(ens_sd);
        rms(time_ind) = rms_err;
   sd_final(time_ind) = sd_mean;

end % End of time for loop


%----------------------------------------------------------------------
% recover some global attributes for plotting purposes.
%----------------------------------------------------------------------

truefile    = netcdf(fname,'nowrite');   % open with low-level netcdf operators.
assimstring = truefile.assim_model_source(:);   % grab a global attribute
modelstring = truefile.model_source(:);         % ditto

modelpath = sscanf(modelstring,'%*s %s %*s');
ind = find(modelpath == filesep);

if isempty(ind)
   disp('do nothing')
else
   for i = 1:length(ind)-1,
      % just the nodes, do not include file separators.
      substr = modelpath(ind(i)+1:ind(i+1)-1);
      switch lower(substr)
         case 'models'   % then the actual model is the next one.
         model = modelpath(ind(i+1)+1:ind(i+2)-1);
      otherwise
         % do nothing
      end
   end
end

disp(sprintf('Selected model is %s',model))

%Open a file for output of summary statistics from second half of run

fid = fopen(fullfile(dir_name,'summary.ascii'), 'w');

fprintf(fid, ['Summary data for  ', dir_name, '\n']);

% Output some graphics of error as function of lead
figure(1);
   subplot(2, 1, 1);
   hold on;
   plot(rms(:));
   plot(sd_final(:), '--');
   title(sprintf('%s %s %s',model, dir_name,': Domain Average RMS Error']);
   grid on;

subplot(2, 1, 2);
hold on;
plot(rms(num_times / 2 : num_times));
plot(sd_final(num_times / 2 : num_times), '--');
grid on;
% Compute the means over this time interval
mean_rms = mean(rms(num_times / 2 : num_times))
mean_sd  = mean(sd_final(num_times / 2: num_times))

legend( strcat('mean rms = ', num2str(mean_rms)), ...
        strcat('mean sd  = ', num2str(mean_sd )),0)

fprintf(fid, ['P ', num2str(mean_rms), '\n']);


%
%
%





