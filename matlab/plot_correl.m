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
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist('diagn_file') ~=1)
   disp('Input name of prior or posterior diagnostics file;')
   diagn_file = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(diagn_file)
      diagn_file = 'Prior_Diag.nc';
   end
end 

vars = CheckModel(diagn_file);
pinfo.fname = diagn_file;

switch lower(vars.model)
   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04', ...
	 'forced_lorenz_96','ikeda'}

      pinfo.base_var = vars.def_var;

      inputstring = input( ...
           sprintf('Input index for base variable (between %d and %d)  ', ...
           vars.min_state_var,vars.max_state_var),'s');
      pinfo.base_var_index = str2num(deblank(inputstring));

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d)  ', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      disp(sprintf('Using diagnostic file %s',diagn_file))
      disp(sprintf('Correlating variable %s index %d at time %d.', ...
           pinfo.base_var, pinfo.base_var_index, pinfo.base_time))

   case {'lorenz_96_2scale'}

      disp(sprintf('Your choice of variables is ''X'' or ''Y'''))
      disp(sprintf('''X'' can range from %d to %d', vars.min_X_var, vars.max_X_var))
      disp(sprintf('''Y'' can range from %d to %d', vars.min_Y_var, vars.max_Y_var))

      % parsing the result of this one is a bit tricky.
      inputstring = input('Input base variable and index i.e.  X 5\n','s');
      [pinfo.base_var, pinfo.base_var_index] = ParseAlphaNumerics(inputstring);

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d)  ', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      disp(sprintf('Using diagnostic file %s',diagn_file))
      disp(sprintf('Correlating variable %s index %d at time %d.', ...
           pinfo.base_var,pinfo.base_var_index, pinfo.base_time))

   case {'simple_advection'}

      disp('Your choice of variables are:')
      disp(vars.vars)
      disp(sprintf('the indices (locations) can range from %d to %d', ...
           vars.min_state_var, vars.max_state_var))

      str1 = sprintf('Input base variable and index i.e. %s %d\n', ...
                      vars.def_var,vars.def_state_vars(1));
      inputstring = input(str1,'s');
      [pinfo.base_var, pinfo.base_var_index] = ParseAlphaNumerics(inputstring);

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d)  ', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      disp(sprintf('Using diagnostic file %s',diagn_file))
      disp(sprintf('Correlating variable %s index %d at time %d.', ...
           pinfo.base_var,pinfo.base_var_index, pinfo.base_time))

   case 'fms_bgrid'

      pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotCorrel');

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

pinfo

PlotCorrel( pinfo );
clear vars inputstring inds str1 vrbl vrbl_inds
