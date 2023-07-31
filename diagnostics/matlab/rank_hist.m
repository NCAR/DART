function bins = rank_hist(ens, verif)
%% RANK_HIST: private function to compute a rank histogram given time series of ensemble and verification
%
% The function to create and plot the rank histogram is called "plot_bins".
%
% Example
%
% truth_file = 'true_state.nc';
% diagn_file = 'preassim.nc';
% plot_bins

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

[ens_size, num_times] = size(ens);

% If one of the dimensions of 'ens' is unity - it must be time.
% If that happens, transpose these things. Why can't Matlab just
% leave the dimensions alone! There has already been a check to make
% sure the ensemble size is greater than 2, so it cannot be 1 here.

if (ens_size == 1)
    ens = ens';
    num_times = size(ens,1);
    ens_size  = size(ens,2);
end

% Need ens_size + 1 bins
bins(1:ens_size + 1) = 0.0;

% Doing this for possibly duplicate ensemble members
for time = 1:num_times
   sens = sort(ens(:, time));
   v = verif(time);
   if(v < sens(1)) 
      bins(1) = bins(1) + 1;
   elseif(v > sens(ens_size))
      bins(ens_size + 1) = bins(ens_size + 1) + 1;
   end
   for i = 1:ens_size - 1
      % Is it between two of the sorted members; this is the traditional easy case
      if(v > sens(i) && v < sens(i + 1))
         bins(i + 1) = bins(i + 1) + 1;
      end
   end

   % What about the case of exactly equal and the possibility of duplicates?
   exact_count = 0;
   for i = 1:ens_size
      if(v == sens(i)) 
         exact_count = exact_count + 1;
      end
   end
   % Each bin that has an endpoint equal to the ensemble member gets 1 / (exact_count + 1) mass
   is_first = true;
   for i = 1:ens_size
      if(v == sens(i))
         if(is_first) 
            bins(i) = bins(i) + 1 / (exact_count + 1);
            is_first = false;
         end 
         bins(i+1) = bins(i+1) + 1 / (exact_count + 1);
      end
   end

end



return

% Loop through time series to get count for each bin
for itime = 1:num_times
   count = 0;
   for imem = 1:ens_size
      if verif(itime) > ens(imem,itime)
         count = count + 1;
      end
   end
   bins(count+1) = bins(count+1) + 1;
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
