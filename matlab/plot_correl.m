% DART : Plots space-time series of correlation between a given variable 
%               at a given time and other variables at all times in an 
%               ensemble time sequence.
%
% plot_correl  interactively queries for the information needed to create
%              the desired correlations.
%              Since different models potentially need different pieces 
%              of information ... the model types are determined and 
%              additional user input may be queried.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end                                                                          
end 

vars = CheckModel(diagn_file);   % also gets default values for this file.

switch lower(vars.model)
   case {'9var','lorenz_63','lorenz_96'}

      pinfo.fname = diagn_file;

      inputstring = input( ...
           sprintf('Input index for base variable (between %d and %d)  ', ...
           vars.min_state_var,vars.max_state_var),'s');
      pinfo.base_var_index = str2num(deblank(inputstring));

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d)  ', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      disp(sprintf('Using diagnostic file %s',diagn_file))
      disp(sprintf('Correlating state variable %d at time %d.', ...
           pinfo.base_var_index, pinfo.base_time))

   case 'fms_bgrid'

      pinfo = GetBgridInfo(diagn_file, 'PlotCorrel');
      pinfo.fname = diagn_file;

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

PlotCorrel( pinfo );
clear vars inputstring
