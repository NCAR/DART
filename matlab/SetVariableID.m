function varid = SetVariableID(vars);
% SetVariableID   queries the 
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

switch lower(vars.model)

   case 'fms_bgrid'

      varid = [1 2 3 4];

   case 'lorenz_96'

      % query to see if these are OK, if not ...

      disp(sprintf('\nUsing state variable IDs %s ', ...
                      num2str(vars.def_state_vars)))
      disp('If these are OK, <cr>;')
      disp(sprintf('If not, please enter array of state variable ID''s (between %d and %d)',vars.min_state_var, vars.max_state_var))
      IDstring = input('(no syntax required)\n','s');
      if isempty(IDstring)
         varid = vars.def_state_vars;
      else
         varid = str2num(IDstring);
      end 

   otherwise

      % low-order models have no choice.
      varid = vars.def_state_vars;

end
