function varid = SetVariableID(vars);
%% SetVariableID   queries the user to override default variable ID's 
%                 (i.e. model state variable indices) for different model types. 
%
% varid = SetVariableID(vars);
%
% This routine is not intended to be called directly.
% As such, the example seems pretty silly.
%
% vars.model = 'lorenz_96';
% vars.def_var = 'state';
% vars.def_state_vars = [2 4 5];
% varid = SetVariableID(vars) 

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

switch lower(vars.model)

   case 'fms_bgrid'

      varid = [1 2 3 4];

   case 'pe2lyr'

      %  model has no choice.
      varid = vars.vars;

   case 'forced_lorenz_96'

      % Determine which State Variable and where ...

      fprintf('\nUsing %s state variable IDs %s\n', ...
                      vars.def_var, num2str(vars.def_state_vars))
      disp('If these are OK, <cr>;')
      disp('If not, please enter array of state variable ID''s')
      fprintf('To choose from entire state enter A 25 50 75 (between %d and %d)\n',vars.min_state_var,vars.max_state_var)
      fprintf('To choose traditional model state enter S 1  23 40 (between %d and %d)\n',vars.min_model_var,vars.max_model_var)
      fprintf('To choose forcing estimates enter F 2 12 22 (between %d and %d)\n',vars.min_force_var,vars.max_force_var)
      IDstring = input('(no intervening syntax required)\n','s');

      if isempty(IDstring) 
         varid = struct('var',vars.def_var,'var_inds',vars.def_state_vars); 
      else 
         [vrbl, vrbl_inds] = ParseAlphaNumerics(IDstring);

	 % Must shift the 'forcing' variables by the number of state variables.
         switch lower(vrbl)
	 case 'f'
            varid = struct('var','state','var_inds',vrbl_inds+vars.num_model_vars); 
	 otherwise
            varid = struct('var','state','var_inds',vrbl_inds); 
	 end
      end 

   case 'lorenz_96_2scale'

      % query to see if these are OK, if not ...

      fprintf('\nUsing %s state variable IDs %s\n', ...
                      vars.def_var, num2str(vars.def_state_vars))
      disp('If these are OK, <cr>;')
      disp('If not, please enter array of state variable ID''s')
      fprintf('To choose X (slow) variables enter X 1  23 40 (between %d and %d)\n',vars.min_X_var,vars.max_X_var)
      fprintf('To choose Y (fast) variables enter Y 23 34 89 (between %d and %d)\n',vars.min_Y_var,vars.max_Y_var)
      IDstring = input('(no intervening syntax required)\n','s');

      if isempty(IDstring) 
         varid = struct('var',vars.def_var,'var_inds',vars.def_state_vars); 
      else 
         [vrbl, vrbl_inds] = ParseAlphaNumerics(IDstring);
         varid = struct('var',upper(vrbl),'var_inds',vrbl_inds); 
      end 

   case 'simple_advection'

      % query to see if the defaults are OK ...

      fprintf('\nUsing variable ''%s'' IDs %s\n', ...
                   vars.def_var, num2str(vars.def_state_vars))
      disp('If these are OK, <cr>;')

      disp('Possible choices of variables are ')
      disp(vars.vars)

      fprintf('and a range between %d and %d\n', vars.min_state_var,  ...
                                                    vars.max_state_var)
      fprintf('you can choose something like ''%s %d %d'', for example.\n', ...
                   vars.vars{1},vars.def_state_vars(1),vars.def_state_vars(2))

      IDstring = input('(no intervening syntax required)\n','s');

      if isempty(IDstring) 
         varid = struct('var',vars.def_var,'var_inds',vars.def_state_vars); 
      else 
         [vrbl, vrbl_inds] = ParseAlphaNumerics(IDstring);
         varid = struct('var',vrbl,'var_inds',vrbl_inds); 
      end 

   case {'lorenz_96','lorenz_04'}

      % query to see if these are OK, if not ...

      fprintf('\nUsing state variable IDs %s\n ', ...
                      num2str(vars.def_state_vars))
      disp('If these are OK, <cr>;')
      fprintf('If not, please enter array of state variable ID''s (between %d and %d)\n',vars.min_state_var, vars.max_state_var)
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
