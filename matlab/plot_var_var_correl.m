% plot_var_var_correl.m
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.

if (exist('fname') ~=1)
   disp('Input name of file;')
   fname = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(fname)
      fname = 'Prior_Diag.nc';
   end                                                                          
end 

vars = CheckModel(fname);   % also gets default values for this file.

inputstring = input(sprintf('Input index for base variable (between %d and %d)  ', ...
                               vars.min_state_var, vars.max_state_var), 's');
base_var_index = str2num(deblank(inputstring));

inputstring = input(sprintf('Input time index for base point (between 1 and %d)  ', ...
                               vars.time_series_length),'s');
base_time = str2num(deblank(inputstring));

inputstring = input(sprintf('Input variable index for correlation (between %d and %d)  ', ...
                               vars.min_state_var, vars.max_state_var), 's');
state_var_index = str2num(deblank(inputstring));

disp(sprintf('Using diagnostic file %s',diagn_file))
disp(sprintf('Correlating state variable %d at time %d with state variable %d.', ...
              base_var_index, base_time, state_var_index))


% could/should check input for valid range, etc.

PlotVarVarCorrel(diagn_file, base_var_index, base_time, state_var_index)
clear vars inputstring
