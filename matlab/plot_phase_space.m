% plot_phase_space.m
% Plots a 3D trajectory of (3 state variables of) a single ensemble member.
% It is possible to overlay subsequent trajectories as follows:
%
% fname = 'True_State.nc'
% var1  = 1;         % variable ID to be used as 'X'
% var2  = 2;         % variable ID to be used as 'Y'
% var3  = 3;         % variable ID to be used as 'Z'
% ens_mem_id = 1;       % ensemble member ID 
% ltype = 'b-';      % line type ('help plot' for details)
% plot_phase_space
%
% hold on
% fname = 'Posterior_Diag.nc'
% var1  = 1;         % variable ID to be used as 'X'
% var2  = 2;         % variable ID to be used as 'Y'
% var3  = 3;         % variable ID to be used as 'Z'
% ens_mem_id = 7;       % ensemble member ID 
% ltype = 'r-';      % line type ('help plot' for details)
% plot_phase_space


if (exist('fname') ~=1)
   disp('Input name of netCDF file;')
   fname = input('<cr> for True_State.nc\n','s');
   if isempty(fname)
      fname = 'True_State.nc';
   end                                                                          
end 


if (exist('var1') ~=1)
   disp('Input state variable index for ''X'' variable.')
   inputstring = input('<cr> for 1\n','s');
   if isempty(inputstring)
      var1 = 1;
   else
      var1 = str2num(deblank(inputstring));
   end                                                                          
end 


if (exist('var2') ~=1)
   disp('Input state variable index for ''Y'' variable.')
   inputstring = input('<cr> for 2\n','s');
   if isempty(inputstring)
      var2 = 2;
   else
      var2 = str2num(deblank(inputstring));
   end                                                                          
end 


if (exist('var3') ~=1)
   disp('Input state variable index for ''Y'' variable.')
   inputstring = input('<cr> for 3\n','s');
   if isempty(inputstring)
      var3 = 3;
   else
      var3 = str2num(deblank(inputstring));
   end                                                                          
end 


if (exist('ens_mem_id') ~=1)
   disp('Input ensemble member ID.')
   inputstring = input('<cr> for 1\n','s');
   if isempty(inputstring)
      ens_mem_id = 1;
   else
      ens_mem_id = str2num(deblank(inputstring));
   end                                                                          
end 


if (exist('ltype') ~=1)
   disp('Input line type string.')
   inputstring = input('<cr> for ''k-''\n','s');
   if isempty(fname)
      ltype = 'k-';
   else
      ltype = inputstring;
   end                                                                          
end 

disp(sprintf('Using file %s, ensemble member %d.',fname,ens_mem_id))
disp(sprintf('Plotting state variables %d %d %d with line type %s.', ...
              var1, var2, var3, ltype))

h = PlotPhaseSpace(fname, ens_mem_id, var1, var2, var3, ltype);
