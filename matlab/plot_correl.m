% plot_correl   Plots space-time series of correlation between a given variable at a given
% time and all other variable at all times in an ensemble time sequence

if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end                                                                          
end 

vars = CheckModel(diagn_file);   % also gets default values for this file.

inputstring = input(sprintf('Input index for base variable (between %d and %d)  ', ...
                               vars.min_state_var,vars.max_state_var),'s');
base_var_index = str2num(deblank(inputstring));

inputstring = input(sprintf('Input time index for base point (between 1 and %d)  ', ...
                               vars.time_series_length),'s');
base_time = str2num(deblank(inputstring));

disp(sprintf('Using diagnostic file %s',diagn_file))
disp(sprintf('Correlating state variable %d at time %d.', ...
              base_var_index, base_time))

PlotCorrel(diagn_file, base_var_index, base_time)

clear vars inputstring
