function PlotVarVarCorrel( pinfo )
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.
% The correlation is done across ensemble members.
%
% PlotVarVarCorrel is intended to be called by 'plot_var_var_correl'
%
% USAGE: PlotVarVarCorrel(fname, base_var_index, base_time, state_var_index)
%
% fname
% base_var_index    index of one of the state variables to correlate  
% base_time         index of the time of interest (timestep)
% state_var_index   index of the other state variable of interest.
%
% Example  (lorenz 63 model with 1000 timesteps)
%%--------------------------------------------------------
% pinfo.fname = 'Posterior_Diag.nc';
% pinfo.base_var_index    = 2;
% pinfo.base_time         = 500;
% pinfo.state_var_index   = 1;
% PlotVarVarCorrel( pinfo )

% TJH Wed Jul  2 09:52:18 MDT 2003

if (exist(pinfo.fname) ~= 2), error(sprintf('%s does not exist.',pinfo.fname)), end

% disp(sprintf('PlotVarVarCorrel: fname is %s',pinfo.fname))
% disp(sprintf('PlotVarVarCorrel: base_var_index is %d',pinfo.base_var_index))
% disp(sprintf('PlotVarVarCorrel: base_time      is %d',pinfo.base_time))

% Get some file-specific information.
f = netcdf(pinfo.fname,'nowrite');
model      = f.model(:);
num_vars   = ncsize(f('StateVariable')); % determine # of state variables
num_times  = ncsize(f('time')); % determine # of output times
num_copies = ncsize(f('copy')); % determine # of ensemble members
close(f);

disp(sprintf('PlotVarVarCorrel: num_vars is %d',num_vars))

% The Base Variable Index must be a valid state variable

if ( pinfo.base_var_index > num_vars )
   disp( sprintf('%s only has %d state variables', pinfo.fname, num_vars))
   error(sprintf('you wanted variable # %d ', pinfo.base_var_index))
end

% The Time must be within range also.

if ( pinfo.base_time > num_times )
   disp( sprintf('%s only has %d output times', pinfo.fname, num_times))
   error(sprintf('you wanted time # %d ', pinfo.base_time))
end

% The State Variable Index must also be a valid state variable

if ( pinfo.base_var_index > num_vars )
   disp( sprintf('%s only has %d state variables', pinfo.fname, num_vars))
   error(sprintf('you wanted variable # %d ', pinfo.base_var_index))
end

% Get 'standard' ensemble series 
 base_var = get_ens_series(pinfo.fname,  pinfo.base_var_index);
state_var = get_ens_series(pinfo.fname, pinfo.state_var_index);

% perform a single correlation
correl = ens_correl(base_var, pinfo.base_time, state_var);

clf; plot(correl);

s1 = sprintf('%s Correlation of state variable %d, T = %d, with variable %d', ...
         model, pinfo.base_var_index, pinfo.base_time, pinfo.state_var_index);
s2 = sprintf('%d ensemble members -- %s', num_copies-2,pinfo.fname); 
title({s1,s2},'interpreter','none','fontweight','bold')
xlabel('time (timestep #)')
ylabel('correlation')

% call out the time index in question, and put a corr==0 reference line.
ax = axis;
hold on;
plot([pinfo.base_time pinfo.base_time],[ -1 1 ],'k:', ...
     [ax(1)         ax(2)],[  0 0 ],'k:')

%axis(ax)
