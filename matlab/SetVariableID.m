function varid = SetVariableID(vars);
% SetVariableID   queries the 

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

switch lower(vars.model)

   case 'fms_bgrid'

      varid = [1 2 3 4];

   case 'pe2lyr'

      %  model has no choice.
      varid = vars.def_state_vars;

   case 'lorenz_96_2scale'

      % query to see if these are OK, if not ...

      disp(sprintf('\nUsing %s state variable IDs %s ', ...
                      vars.def_var, num2str(vars.def_state_vars)))
      disp('If these are OK, <cr>;')
      disp('If not, please enter array of state variable ID''s')
      disp(sprintf('To choose X (slow) variables enter X 1  23 40 (between %d and %d)',vars.min_X_var,vars.max_X_var))
      disp(sprintf('To choose Y (fast) variables enter Y 23 34 89 (between %d and %d)',vars.min_Y_var,vars.max_Y_var))
      IDstring = input('(no intervening syntax required)\n','s');

      if isempty(IDstring) 
         varid = struct('var',vars.def_var,'var_inds',vars.def_state_vars); 
      else 
         [vrbl, vrbl_inds] = ParseAlphaNumeric(IDstring);
         varid = struct('var',vrbl,'var_inds',vrbl_inds); 
      end 

   case {'lorenz_96','lorenz_04'}

      % query to see if these are OK, if not ...

      disp(sprintf('\nUsing state variable IDs %s ', ...
                      num2str(vars.def_state_vars)))
      disp('If these are OK, <cr>;')
      disp(sprintf('If not, please enter array of state variable ID''s (between %d and %d)',vars.min_state_var, vars.max_state_var))
      IDstring = input('(no syntax required)\n','s');
      if isempty(IDstring)
         varid = struct('var',vars.def_var,'var_inds',vars.def_state_vars); 
      else
         varid = struct('var',vars.def_var,'var_inds',str2num(IDstring)); 
      end 

   otherwise

      % ultra low-order models have no choice.
      varid = struct('var',vars.def_var,'var_inds',vars.def_state_vars); 

end
