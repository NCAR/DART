function PlotCorrel(fname, base_var_index, base_time)
% PlotCorrel   space-time series of correlation between a variable at a given
% time and all variables at all times in an ensemble time sequence.
%
% PlotCorrel is intended to be called by 'plot_correl'.
%
% USAGE: PlotCorrel(fname, base_var_index, base_time)
%
% fname             name of netCDF file containing a DART ensemble
% base_var_index    index of state variable used as standard in correlation
% base_time         index of time series to use as the standard for correlation
%
% Example 1   (9var model with 1000 time steps)
%%------------------------------------------------------------------
% fname          = 'Prior_Diag.nc';
% base_var_index = 5;                           % picked arbitrarily
% base_time      = 238;                         % ditto
% PlotCorrel(fname, base_var_index, base_time)  % generates a plot

% TJH Wed Jul  2 08:39:46 MDT 2003

if (exist(fname) ~= 2), error(sprintf('%s does not exist.',fname)), end

% Get some file-specific information.
f = netcdf(fname,'nowrite');
model      = f.model(:);
num_vars   = ncsize(f{'StateVariable'}); % determine # of state variables
num_times  = ncsize(f{'time'}); % determine # of output times
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
close(f)

% The Base Variable Index must be a valid state variable
if ( base_var_index > num_vars )
   disp( sprintf('%s only has %d state variables', fname, num_vars))
   error(sprintf('you wanted variable # %d ', base_var_index))
end

% The Time must be within range also.
if ( base_time > num_times )
   disp( sprintf('%s only has %d output times', fname, num_times))
   error(sprintf('you wanted time # %d ', base_time))
end

statevariables = getnc(fname,'StateVariable');

% Get 'standard' ensemble series 
base = get_ens_series(fname, base_var_index);

% It is efficient to preallocate correl storage ... 
correl = zeros(num_vars,num_times);

% Need to loop through all variables in the ensemble
for i = 1:num_vars
   state_var = get_ens_series(fname, i);
   correl(i, :) = ens_correl(base, base_time, state_var);
end

% Now for the plotting part ...
clf;

contour(correl,[-1:0.2:1]);
s1 = sprintf('%s Correlation of state variable %d, T = %d of %s', ...
         model, base_var_index, base_time,fname);
s2 = sprintf('against all variables, all times, all %d ensemble members', ...
         num_copies-2); 
title({s1,s2},'interpreter','none','fontweight','bold')
xlabel('time (timestep #)')
ylabel('state variable (index)')
set(gca,'YTick',statevariables)
colorbar

% highlight the reference state variable and time

hold on;
plot(base_time,base_var_index,'kh','MarkerSize',12,'MarkerFaceColor','k')
