% plot_phase_space.m
% Plots a 3D trajectory of (3 state variables of) a single ensemble member.
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
% Because it IS possible to overlay plots, the onus is on YOU to make
% sure the current figure is "cleared" before you plot the first
% trajectory.

% TJH Wed Jul  2 10:11:04 MDT 2003


if (exist('fname') ~=1)
   fname = input('Input name of netCDF file; <cr> for True_State.nc  ','s');
   if isempty(fname)
      fname = 'True_State.nc';
   end                                                                          
end 
pinfo.fname = fname;

vars  = CheckModel(fname);   % also gets default values for this model.
varid = SetVariableID(vars);      % queries for variable IDs if needed.

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_96'}

      if (isfield(pinfo,'var1') ~=1)
         s1 = input('Input state variable index for ''X'' variable. <cr> for 1  ','s');
         if isempty(s1), pinfo.var1 = 1; else pinfo.var1 = str2num(deblank(s1)); end
      end 

      if (isfield(pinfo,'var2') ~=1)
         s1 = input('Input state variable index for ''Y'' variable. <cr> for 2  ','s');
         if isempty(s1), pinfo.var2 = 2; else pinfo.var2 = str2num(deblank(s1)); end
      end 

      if (isfield(pinfo,'var3') ~=1)
         s1 = input('Input state variable index for ''Z'' variable. <cr> for 3  ','s');
         if isempty(s1), pinfo.var3 = 3; else pinfo.var3 = str2num(deblank(s1)); end
      end 

      if (isfield(pinfo,'ens_mem') ~=1)
         s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
         if isempty(s1), pinfo.ens_mem = 'true state'; else pinfo.ens_mem = s1; end
      end 

      if (isfield(pinfo,'ltype') ~=1)
         s1 = input('Input line type string. <cr> for ''k-''  ','s');
         if isempty(fname), pinfo.ltype = 'k-'; else pinfo.ltype = s1; end
      end 
      
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
