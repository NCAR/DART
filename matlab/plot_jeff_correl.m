% plot_jeff_correl.m
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.

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

if (exist('fname') ~=1)
   disp('Input name of file;')
   fname = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(fname)
      fname = 'Prior_Diag.nc';
   end
else
   % check to make sure they are using a file with some ensemble members.

   f          = netcdf(fname,'nowrite');
   var_atts   = dim(f{'copy'});
   num_copies = length(var_atts{1}); % determine # of ensemble members
   close(f)

   if (num_copies <= 3) 
      error(sprintf('Sorry -- %s does not have enough ensemble members to correlate.',fname))
   end 
end 

vars = CheckModel(fname);   % also gets default values for this file.
pinfo.fname = fname;

switch lower(vars.model) 
   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04', ...
	 'forced_lorenz_96','ikeda'}

      pinfo.base_var  = vars.def_var;
      pinfo.state_var = vars.def_var;

      inputstring = input( ...
           sprintf('Input index for base variable (between %d and %d, inclusive)  ', ...
           vars.min_state_var, vars.max_state_var), 's');
      pinfo.base_var_index = str2num(deblank(inputstring));

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d, inclusive)  ', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      inputstring = input( ...
           sprintf('Input variable index for correlation (between %d and %d, inclusive)  ', ...
           vars.min_state_var, vars.max_state_var), 's');
      pinfo.state_var_index = str2num(deblank(inputstring));

  %   disp(sprintf('Using diagnostic file %s',fname))
  %   disp(sprintf('Correlating variable %s %d at time %d with variable %s %d.', ...
  %           pinfo.def_var, pinfo.base_var_index, pinfo.base_time, pinfo.state_var, pinfo.state_var_index))

   case {'lorenz_96_2scale'}

      disp(sprintf('Your choice of variables is ''X'' or ''Y'''))
      disp(sprintf('''X'' can range from %d to %d', vars.min_X_var, vars.max_X_var))
      disp(sprintf('''Y'' can range from %d to %d', vars.min_Y_var, vars.max_Y_var))

      inputstring = input('Input base variable and index i.e.  X 5\n','s');
      [pinfo.base_var, pinfo.base_var_index] = ParseAlphaNumerics(inputstring);

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d, inclusive) \n', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      inputstring = input('Input variable and index for correlation \n','s');
      [pinfo.state_var, pinfo.state_var_index] = ParseAlphaNumerics(inputstring);

   case {'simple_advection'}

      disp('Your choice of variables are:')
      disp(vars.vars)
      disp(sprintf('the indices (locations) can range from %d to %d, inclusive', ...
           vars.min_state_var, vars.max_state_var))

      str1 = sprintf('Input base variable and index i.e. %s %d\n', ...
                      vars.def_var,vars.def_state_vars(1));
      inputstring = input(str1,'s');
      [pinfo.base_var, pinfo.base_var_index] = ParseAlphaNumerics(inputstring);

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d, inclusive)  ', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      inputstring = input('Input variable and index for correlation \n','s');
      [pinfo.state_var, pinfo.state_var_index] = ParseAlphaNumerics(inputstring);

   case 'fms_bgrid'

      pinfo = GetBgridInfo(pinfo, fname, 'PlotVarVarCorrel');

   case 'pe2lyr'

      pinfo = GetPe2lyrInfo(pinfo, fname, 'PlotVarVarCorrel');

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

% could/should check input for valid range, etc.

PlotJeffCorrel( pinfo )
clear vars inputstring fname
