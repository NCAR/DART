% plot_phase_space : Plots a 3D trajectory of (3 state variables of) a single ensemble member.
%
% Because it IS possible to overlay plots, the onus is on YOU to make
% sure the current figure is "cleared" before you plot the first
% trajectory.
%
% It is possible to overlay subsequent trajectories as follows:
%
% clf;                      % clears the current figure  
% fname = 'True_State.nc';
% var1  = 1;                % variable ID to be used as 'X'
% var2  = 2;                % variable ID to be used as 'Y'
% var3  = 3;                % variable ID to be used as 'Z'
% ens_mem = 'true state';   % ensemble member metadata string
% ltype = 'b-';             % line type ('help plot' for details)
% plot_phase_space
%
% hold on
% fname      = 'Posterior_Diag.nc';
% ens_mem    = 'ensemble mean';           % ensemble member ID 
% ltype      = 'r-';        % line type ('help plot' for details)
% plot_phase_space
%
% ens_mem    = 'ensemble member4';        % ensemble member ID 
% ltype      = 'c-';        % line type ('help plot' for details)
% plot_phase_space
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist('fname','var') ~=1)
   fname = input('Input name of netCDF file; <cr> for True_State.nc  ','s');
   if isempty(fname)
      fname = 'True_State.nc';
   end                                                                          
else
   if isempty(fname), fname = 'True_State.nc'; end
   s1 = input(sprintf('Input name of netCDF file. <cr> for  %s ',fname),'s');
   if ~isempty(s1), fname = s1; end
end 

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

pinfo.fname = strtrim(fname);

vars  = CheckModel(fname);   % also gets default values for this model.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04','forced_lorenz_96'}

      if (ishold), clear var1 var2 var3 ens_mem ltype; end

      str1 = sprintf('[%d - %d]',vars.min_state_var, vars.max_state_var);

      if (exist('var1') ~=1)
         s1 = input(sprintf('Input variable index for ''X'' variable %s. <cr> for 1.  ',str1),'s');
         if isempty(s1), var1 = 1; else var1 = str2num(deblank(s1)); end
      end 

      if (exist('var2') ~=1)
         s1 = input(sprintf('Input variable index for ''Y'' variable %s. <cr> for 2.  ',str1),'s');
         if isempty(s1), var2 = 2; else var2 = str2num(deblank(s1)); end
      end 

      if (exist('var3') ~=1)
         s1 = input(sprintf('Input variable index for ''Z'' variable %s. <cr> for 3.  ',str1),'s');
         if isempty(s1), var3 = 3; else var3 = str2num(deblank(s1)); end
      end 

      if (exist('ens_mem') ~=1)
         s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
         if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end
      end 

      if (exist('ltype') ~=1)
         s1 = input('Input line type string. <cr> for ''k-''  ','s');
         if isempty(s1), ltype = 'k-'; else ltype = s1; end
      end 

      pinfo = struct('fname'   , fname       , ...
                     'model'   , vars.model  , ...
                     'var1name', vars.def_var, 'var1ind', var1, ...
                     'var2name', vars.def_var, 'var2ind', var2, ...
                     'var3name', vars.def_var, 'var3ind', var3, ...
                     'ens_mem' , ens_mem     , ...
                     'ltype'   , ltype   );

   case {'lorenz_96_2scale'}

      if (ishold), clear var1 var2 var3 ens_mem ltype; end

      fprintf('Your choice of variables is ''X'' or ''Y''\n')
      fprintf('''X'' can range from %d to %d\n', vars.min_X_var, vars.max_X_var)
      fprintf('''Y'' can range from %d to %d\n', vars.min_Y_var, vars.max_Y_var)

      % really should utilize the defaults ... but its getting late.
      inputstring = input('Input variable and index for axis 1 i.e.  X 5\n','s');
      [var1name, var1] = ParseAlphaNumerics(inputstring);

      inputstring = input('Input variable and index for axis 2 i.e.  X 18\n','s');
      [var2name, var2] = ParseAlphaNumerics(inputstring);

      inputstring = input('Input variable and index for axis 3 i.e.  X 24\n','s');
      [var3name, var3] = ParseAlphaNumerics(inputstring);

      if (exist('ens_mem') ~=1)
         disp('It is necessary to pick an ensemble member to plot.')
         disp('Since we pick it based on the metadata string, it could be:')
         disp('''true_state'', ''ensemble mean'', ''ensemble member10'' ... you get it.')
         s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
         if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end
      end 

      if (exist('ltype') ~=1)
         s1 = input('Input line type string. <cr> for ''k-''  ','s');
         if isempty(s1), ltype = 'k-'; else ltype = s1; end
      end 

      pinfo = struct('fname'   , fname    , ...
                     'model'   , vars.model  , ...
                     'var1name', var1name , 'var1ind', var1, ...
                     'var2name', var2name , 'var2ind', var2, ...
                     'var3name', var3name , 'var3ind', var3, ...
                     'ens_mem' , ens_mem  , ...
                     'ltype'   , ltype   );

  case {'simple_advection'}

      if (ishold), clear var1 var2 var3 ens_mem ltype; end

      disp('Your choice of variables are:')
      disp(vars.vars)
      fprintf('the indices (locations) can range from %d to %d\n', ...
           vars.min_state_var, vars.max_state_var)

      str1 = sprintf('Input variable and index for axis 1 <cr> for %s %d\n', ...
                      vars.def_var,vars.def_state_vars(1));
      str2 = sprintf('Input variable and index for axis 2 <cr> for %s %d\n', ...
                      vars.def_var,vars.def_state_vars(2));
      str3 = sprintf('Input variable and index for axis 2 <cr> for %s %d\n', ...
                      vars.def_var,vars.def_state_vars(3));

      s1 = input(str1,'s');
      if isempty(s1)
         var1name = vars.def_var;
         var1     = vars.def_state_vars(1);
      else
         [var1name, var1] = ParseAlphaNumerics(s1);
      end

      s1 = input(str2,'s');
      if isempty(s1)
         var2name = vars.def_var;
         var2     = vars.def_state_vars(2);
      else
         [var2name, var2] = ParseAlphaNumerics(s1);
      end

      s1 = input(str3,'s');
      if isempty(s1)
         var3name = vars.def_var;
         var3     = vars.def_state_vars(3);
      else
         [var3name, var3] = ParseAlphaNumerics(s1);
      end

      if (exist('ens_mem') ~=1)
         disp('It is necessary to pick an ensemble member to plot.')
         disp('Since we pick it based on the metadata string, it could be:')
         disp('''true_state'', ''ensemble mean'', ''ensemble member10'' ... you get it.')
         s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
         if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end
      end 

      if (exist('ltype') ~=1)
         s1 = input('Input line type string. <cr> for ''k-''  ','s');
         if isempty(s1), ltype = 'k-'; else ltype = s1; end
      end 

      pinfo = struct('fname'   , fname    , ...
                     'model'   , vars.model  , ...
                     'var1name', var1name , 'var1ind', var1, ...
                     'var2name', var2name , 'var2ind', var2, ...
                     'var3name', var3name , 'var3ind', var3, ...
                     'ens_mem' , ens_mem  , ...
                     'ltype'   , ltype   );

   case 'fms_bgrid'

      pinfo = GetBgridInfo(pinfo, fname, 'PlotPhaseSpace');

   case 'cam'

      pinfo = GetCamInfo(pinfo, fname, 'PlotPhaseSpace');

   case 'pe2lyr'

      pinfo = GetPe2lyrInfo(pinfo, fname, 'PlotPhaseSpace');

   case 'mitgcm_ocean'

      pinfo = GetMITgcm_oceanInfo(pinfo, fname, 'PlotPhaseSpace');

   case {'ikeda'}

      if (ishold), clear var1 var2 var3 ens_mem ltype; end

      str1 = sprintf('[%d - %d]',vars.min_state_var, vars.max_state_var);

      if (exist('var1') ~=1)
         s1 = input(sprintf('Input variable index for ''X'' variable %s. <cr> for 1.  ',str1),'s');
         if isempty(s1), var1 = 1; else var1 = str2num(deblank(s1)); end
      end 

      if (exist('var2') ~=1)
         s1 = input(sprintf('Input variable index for ''Y'' variable %s. <cr> for 2.  ',str1),'s');
         if isempty(s1), var2 = 2; else var2 = str2num(deblank(s1)); end
      end 

      if (exist('ens_mem') ~=1)
         s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
         if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end
      end 

      if (exist('ltype') ~=1)
         s1 = input('Input line type string. <cr> for ''k-''  ','s');
         if isempty(s1), ltype = 'k-'; else ltype = s1; end
      end 

      pinfo = struct('fname'   , fname       , ...
                     'model'   , vars.model  , ...
                     'var1name', vars.def_var, 'var1ind', var1, ...
                     'var2name', vars.def_var, 'var2ind', var2, ...
                     'ens_mem' , ens_mem     , ...
                     'ltype'   , ltype   );

   otherwise

      error('model %s not implemented yet', vars.model)

end

PlotPhaseSpace( pinfo );
clear s1
