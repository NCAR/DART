% plot_var_var_correl.m
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.

if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for ../work/Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = '../work/Prior_Diag.nc';
   end                                                                          
end 

inputstring = input('Input index for base variable\n','s');
base_var_index = str2num(deblank(inputstring));

inputstring = input('Input time index for base point\n','s');
base_time = str2num(deblank(inputstring));

inputstring = input('Input variable index for correlation \n','s');
state_var_index = str2num(deblank(inputstring));

disp(sprintf('Using diagnostic file %s',diagn_file))
disp(sprintf('Correlating state variable %d at time %d with state variable %d.', ...
              base_var_index, base_time, state_var_index))

PlotVarVarCorrel(diagn_file, base_var_index, base_time, state_var_index)
