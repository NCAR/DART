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


if (exist('var1') ~=1)
   inputstring = input('Input state variable index for ''X'' variable. <cr> for 1  ','s');
   if isempty(inputstring)
      var1 = 1;
   else
      var1 = str2num(deblank(inputstring));
   end                                                                          
end 


if (exist('var2') ~=1)
   inputstring = input('Input state variable index for ''Y'' variable. <cr> for 2  ','s');
   if isempty(inputstring)
      var2 = 2;
   else
      var2 = str2num(deblank(inputstring));
   end                                                                          
end 


if (exist('var3') ~=1)
   inputstring = input('Input state variable index for ''Z'' variable. <cr> for 3  ','s');
   if isempty(inputstring)
      var3 = 3;
   else
      var3 = str2num(deblank(inputstring));
   end                                                                          
end 


if (exist('ens_mem') ~=1)
   inputstring = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
   if isempty(inputstring)
      ens_mem = 'true state';
   else
      ens_mem = inputstring;
   end                                                                          
end 


if (exist('ltype') ~=1)
   inputstring = input('Input line type string. <cr> for ''k-''  ','s');
   if isempty(fname)
      ltype = 'k-';
   else
      ltype = inputstring;
   end                                                                          
end 

disp(sprintf('Using file %s, ensemble member %s.',fname,ens_mem))
disp(sprintf('Plotting state variables %d %d %d with line type %s.', ...
              var1, var2, var3, ltype))

PlotPhaseSpace(fname, ens_mem, var1, var2, var3, ltype);
clear inputstring
