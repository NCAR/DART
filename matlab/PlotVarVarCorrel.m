function PlotVarVarCorrel(fname, base_var_index, base_time, state_var_index)
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.

if (exist(fname) ~= 2), error(sprintf('%s does not exist.',fname)), end

% disp(sprintf('PlotVarVarCorrel: fname is %s',fname))
% disp(sprintf('PlotVarVarCorrel: base_var_index is %d',base_var_index))
% disp(sprintf('PlotVarVarCorrel: base_time      is %d',base_time))

% Get some file-specific information.
f = netcdf(fname,'nowrite');
model      = f.model(:);
num_vars   = ncsize(f{'StateVariable'}); % determine # of state variables
num_times  = ncsize(f{'time'}); % determine # of output times
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
close(f);

disp(sprintf('PlotVarVarCorrel: num_vars is %d',num_vars))

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

% The State Variable Index must also be a valid state variable

if ( base_var_index > num_vars )
   disp( sprintf('%s only has %d state variables', fname, num_vars))
   error(sprintf('you wanted variable # %d ', base_var_index))
end

% Get 'standard' ensemble series 
 base_var = get_ens_series(fname,  base_var_index);
state_var = get_ens_series(fname, state_var_index);

% perform a single correlation
correl = ens_correl(base_var, base_time, state_var);

clf; plot(correl);

s1 = sprintf('%s Correlation of state variable %d, T = %d, with variable %d', ...
         model, base_var_index, base_time, state_var_index);
s2 = sprintf('%d ensemble members -- %s', num_copies-2,fname); 
title({s1,s2},'interpreter','none','fontweight','bold')
xlabel('time (timestep #)')
ylabel('correlation')

% call out the time index in question, and put a corr==0 reference line.
ax = axis;
hold on;
plot([base_time base_time],[ -1 1 ],'k:', ...
     [ax(1)         ax(2)],[  0 0 ],'k:')

%axis(ax)
