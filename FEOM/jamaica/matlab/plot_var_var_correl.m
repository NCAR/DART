% DART : Plots time series of correlation between 
% a given variable at a given time and another variable at all times 
% in an ensemble time sequence.

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

      inputstring = input( ...
           sprintf('Input index for base variable (between %d and %d)  ', ...
           vars.min_state_var, vars.max_state_var), 's');
      pinfo.base_var_index = str2num(deblank(inputstring));

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d)  ', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));

      inputstring = input( ...
           sprintf('Input variable index for correlation (between %d and %d)  ', ...
           vars.min_state_var, vars.max_state_var), 's');
      pinfo.state_var_index = str2num(deblank(inputstring));

      pinfo.base_var  = vars.def_var;
      pinfo.state_var = vars.def_var;

    % disp(sprintf('Using diagnostic file %s',fname))
    % disp(sprintf('Correlating variable %s %d at time %d with variable %s %d.', ...
    %         pinfo.base_var, pinfo.base_var_index, pinfo.base_time, ...
    %         pinfo.state_var, pinfo.state_var_index))

   case {'lorenz_96_2scale'}

      disp(sprintf('Your choice of variables is ''X'' or ''Y'''))
      disp(sprintf('''X'' can range from %d to %d', vars.min_X_var, vars.max_X_var))
      disp(sprintf('''Y'' can range from %d to %d', vars.min_Y_var, vars.max_Y_var))

      % parsing the result of this one is a bit tricky.
      inputstring = input('Input base variable and index i.e.  X 5\n','s');
      [pinfo.base_var, pinfo.base_var_index] = ParseAlphaNumeric(inputstring);

      inputstring = input( ...
           sprintf('Input time index for base point (between 1 and %d) \n', ...
           vars.time_series_length),'s');
      pinfo.base_time = str2num(deblank(inputstring));


      % parsing the result of this one is a bit tricky.
      inputstring = input('Input variable and index for correlation \n','s');
      [pinfo.state_var, pinfo.state_var_index] = ParseAlphaNumeric(inputstring);

    % disp(sprintf('Using diagnostic file %s',fname))
    % disp(sprintf('Correlating variable %s %d at time %d with variable %s %d.', ...
    %         pinfo.base_var, pinfo.base_var_index, pinfo.base_time, ...
    %         pinfo.state_var, pinfo.state_var_index))

   case 'fms_bgrid'

      pinfo = GetBgridInfo(fname, 'PlotVarVarCorrel');

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

pinfo

PlotVarVarCorrel( pinfo )
clear vars inputstring
