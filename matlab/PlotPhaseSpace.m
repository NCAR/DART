function PlotPhaseSpace(fname,ens_mem,var1,var2,var3,ltype)
% PlotPhaseSpace: Plots trajectories of 3 variables for any ensemble member
%
% PlotPhaseSpace is intended to be called by 'plot_phase_space'
%
% USAGE: PlotPhaseSpace(fname, ens_mem, var1, var2, var3, ltype)
%
% fname        name of netCDF DART file
% ens_mem      ensemble member metadata string
% var1         ID of state variable to plot as 'X'
% var2         ID of state variable to plot as 'Y'
% var3         ID of state variable to plot as 'Z'
% ltype        line type (see 'help plot' for details)
%
% Example 1   ( 9 variable model )
%%--------------------------------------------------------
% fname      = 'True_State.nc';
% ens_mem    = 'true state';    % true state only has 1 ens mem ...
% var1       = 3;    
% var2       = 6;
% var3       = 7;
% ltype      = 'b-'; % solid blue line
% PlotPhaseSpace(fname,ens_mem,var1,var2,var3,ltype)
%
% that worked so well, lets overlay another (using the same state variables)
%
% hold on;    
% fname      = 'Prior_Diag.nc';
% ens_mem    = 'ensemble member4';               % why not?
% ltype      = 'r:';            % plot it in a red 'dotted' line
% PlotPhaseSpace(fname,ens_mem,var1,var2,var3,ltype)
%
% note the legend has both lines annotated.

% TJH Wed Jul  2 09:28:08 MDT 2003

if ( exist(fname) ~= 2 ), error(sprintf('file %s does not exist.',fname)), end

% Get some information from the file 
f = netcdf(fname);
model      = f.model(:);
num_vars   = ncsize(f{'StateVariable'}); % determine # of state variables
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
num_times  = ncsize(f{'time'}); % determine # of output times
close(f);

% rudimentary bulletproofing
if ( (var1 > num_vars) | (var1 < 1) ) 
   disp(sprintf('\n%s has %d state variables',fname,num_vars))
   disp(sprintf('%d  <= ''var1'' <= %d',1,num_vars))
   error(sprintf('var1 (%d) out of range',var1))
end

if ( (var2 > num_vars) | (var2 < 1) ) 
   disp(sprintf('\n%s has %d state variables',fname,num_vars))
   disp(sprintf('%d  <= ''var2'' <= %d',1,num_vars))
   error(sprintf('var2 (%d) out of range',var2))
end

if ( (var3 > num_vars) | (var3 < 1) ) 
   disp(sprintf('\n%s has %d state variables',fname,num_vars))
   disp(sprintf('%d  <= ''var3'' <= %d',1,num_vars))
   error(sprintf('var3 (%d) out of range',var3))
end

ens_mem_id = get_copy_index(fname, ens_mem);  % will error out if no such ens_mem 

x = get_var_series(fname, ens_mem_id, var1);
y = get_var_series(fname, ens_mem_id, var2);
z = get_var_series(fname, ens_mem_id, var3);

% There is no model-dependent segment ...
% As long as you have three variables, this works for all models.

h = plot3(x,y,z,ltype);

% If there is no legend, we are assuming this is the first plot on this
% axis. We need to add a title, legend, axes labels, etc.
[legh, objh, outh, outm] = legend;
if (isempty(legh))

   % title(sprintf('%s ensemble member %d of %s',model,ens_mem_id,fname), ...
   %    'interpreter','none','fontweight','bold')
   title('The Attractor','fontweight','bold')
   xlabel(sprintf('state variable # %d',var1))
   ylabel(sprintf('state variable # %d',var2))
   zlabel(sprintf('state variable # %d',var3))

   s = sprintf('%d %d %d %s %s %d',var1,var2,var3,model,fname,ens_mem_id);
   h = legend(s,0);
   [legh, objh, outh, outm] = legend;
   set(objh(1),'interpreter','none')

else
   % Must add salient information to the legend.
   % legh     handle to the legend axes
   % objh     handle for the text, lines, and patches in the legend
   % outh     handle for the lines and patches in the plot
   % outm     cell array for the text in the legend
   nlines = length(outm);
   outm{nlines+1} = sprintf('%d %d %d %s %s %d', var1, var2, var3, model,...
                             fname, ens_mem_id);
   [legh, objh, outh, outm] = legend([outh; h],outm,0);

   set(objh(1),'interpreter','none')
end
legend boxoff
