% DART : Plots a 3D trajectory of (3 state variables of) a single ensemble member.
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
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% $Source$
% $Revision$
% $Date$

if (exist('fname') ~=1)
   fname = input('Input name of netCDF file; <cr> for True_State.nc  ','s');
   if isempty(fname)
      fname = 'True_State.nc';
   end                                                                          
else
   s1 = input(sprintf('Input name of netCDF file. <cr> for  %s ',fname),'s');
   if ~isempty(s1), fname = str2num(deblank(s1)); end
end 

pinfo.fname = fname;

vars  = CheckModel(fname);   % also gets default values for this model.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04'}

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

      pinfo = struct('fname'  , fname   , ...
                     'var1'   , var1    , ...
                     'var2'   , var2    , ...
                     'var3'   , var3    , ...
                     'ens_mem', ens_mem , ...
                     'ltype'  , ltype   );

      disp(sprintf('Using file %s, ensemble member %s.',pinfo.fname,pinfo.ens_mem))
      disp(sprintf('Plotting state variables %d %d %d with line type %s.', ...
                    pinfo.var1, pinfo.var2, pinfo.var3, pinfo.ltype))

   case 'fms_bgrid'

      pinfo = GetBgridInfo(fname, 'PlotPhaseSpace');

      pinfo                            % just echo stuff for posterity.

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

PlotPhaseSpace( pinfo );
clear s1
