function PlotPhaseSpace( pinfo )
% PlotPhaseSpace: Plots trajectories of 3 variables for any ensemble member
%
% PlotPhaseSpace is intended to be called by 'plot_phase_space'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotPhaseSpace( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% fname        name of netCDF DART file
% ens_mem      ensemble member metadata string
% var1         ID of state variable to plot as 'X'
% var2         ID of state variable to plot as 'Y'
% var3         ID of state variable to plot as 'Z'
% ltype        line type (see 'help plot' for details)
%
% Example 1   ( 9 variable model )
%%--------------------------------------------------------
% pinfo.fname      = 'True_State.nc';
% pinfo.ens_mem    = 'true state';    % true state only has 1 ens mem ...
% pinfo.var1       = 3;    
% pinfo.var2       = 6;
% pinfo.var3       = 7;
% pinfo.ltype      = 'b-'; % solid blue line
% PlotPhaseSpace( pinfo )
%
% that worked so well, lets overlay another (using the same state variables)
%
% hold on;    
% pinfo.fname      = 'Prior_Diag.nc';
% pinfo.ens_mem    = 'ensemble member4';               % why not?
% pinfo.ltype      = 'r:';            % plot it in a red 'dotted' line
% PlotPhaseSpace( pinfo )
%
% note the legend has both lines annotated.

% TJH Wed Jul  2 09:28:08 MDT 2003

if ( exist(pinfo.fname) ~= 2 ), error(sprintf('file %s does not exist.',pinfo.fname)), end

% Get some information from the file 
f          = netcdf(pinfo.fname);
model      = f.model(:);
num_vars   = ncsize(f('StateVariable')); % determine # of state variables
num_copies = ncsize(f('copy')); % determine # of ensemble members
num_times  = ncsize(f('time')); % determine # of output times
close(f);


switch lower(model)

   case {'9var','lorenz_63','lorenz_96'}

      % rudimentary bulletproofing
      if ( (pinfo.var1 > num_vars) | (pinfo.var1 < 1) ) 
         disp(sprintf('\n%s has %d state variables',pinfo.fname,num_vars))
         disp(sprintf('%d  <= ''var1'' <= %d',1,num_vars))
         error(sprintf('var1 (%d) out of range',pinfo.var1))
      end
      
      if ( (pinfo.var2 > num_vars) | (pinfo.var2 < 1) ) 
         disp(sprintf('\n%s has %d state variables',pinfo.fname,num_vars))
         disp(sprintf('%d  <= ''var2'' <= %d',1,num_vars))
         error(sprintf('var2 (%d) out of range',pinfo.var2))
      end
      
      if ( (pinfo.var3 > num_vars) | (pinfo.var3 < 1) ) 
         disp(sprintf('\n%s has %d state variables',pinfo.fname,num_vars))
         disp(sprintf('%d  <= ''var3'' <= %d',1,num_vars))
         error(sprintf('var3 (%d) out of range',pinfo.var3))
      end
      
      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);  % errors out if no ens_mem 
      
      x = get_var_series(pinfo.fname, ens_mem_id, pinfo.var1);
      y = get_var_series(pinfo.fname, ens_mem_id, pinfo.var2);
      z = get_var_series(pinfo.fname, ens_mem_id, pinfo.var3);
      
      % There is no model-dependent segment ...
      % As long as you have three variables, this works for all models.
      
      h = plot3(x,y,z,pinfo.ltype);
      
      % If there is no legend, we are assuming this is the first plot on this
      % axis. We need to add a title, legend, axes labels, etc.
      [legh, objh, outh, outm] = legend;
      if (isempty(legh))
      
         % title(sprintf('%s ensemble member %d of %s',model,ens_mem_id,pinfo.fname), ...
         %    'interpreter','none','fontweight','bold')
         title('The Attractor','fontweight','bold')
         xlabel(sprintf('state variable # %d',pinfo.var1))
         ylabel(sprintf('state variable # %d',pinfo.var2))
         zlabel(sprintf('state variable # %d',pinfo.var3))
      
         s = sprintf('%d %d %d %s %s %s', pinfo.var1, pinfo.var2, ...
                     pinfo.var3, pinfo.fname, model, pinfo.ens_mem);
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
         outm{nlines+1} = sprintf('%d %d %d %s %s %s', ...
              pinfo.var1, pinfo.var2, pinfo.var3, pinfo.fname, ...
              model, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);
      
         set(objh(1),'interpreter','none')
      end
      legend boxoff

   case 'fms_bgrid'

      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);  % errors out if no ens_mem 
      
      x = Get1Copy(pinfo.fname, ens_mem_id, pinfo.var1_var, ...
                  pinfo.var1_lvlind, pinfo.var1_latind, pinfo.var1_lonind);
      y = Get1Copy(pinfo.fname, ens_mem_id, pinfo.var2_var, ...
                  pinfo.var2_lvlind, pinfo.var2_latind, pinfo.var2_lonind);
      z = Get1Copy(pinfo.fname, ens_mem_id, pinfo.var3_var, ...
                  pinfo.var3_lvlind, pinfo.var3_latind, pinfo.var3_lonind);

      % There is no model-dependent segment ...
      % As long as you have three variables, this works for all models.
      
      h = plot3(x,y,z,pinfo.ltype);
      
      % If there is no legend, we are assuming this is the first plot on this
      % axis. We need to add a title, legend, axes labels, etc.
      [legh, objh, outh, outm] = legend;
      if (isempty(legh))
      
         % title(sprintf('%s ensemble member %d of %s',model,ens_mem_id,pinfo.fname), ...
         %    'interpreter','none','fontweight','bold')
         title('The Attractor','fontweight','bold')
         xlabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var1_var, pinfo.var1_lvl, pinfo.var1_lat, pinfo.var1_lon))
         ylabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var2_var, pinfo.var2_lvl, pinfo.var2_lat, pinfo.var2_lon))
         zlabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var3_var, pinfo.var3_lvl, pinfo.var3_lat, pinfo.var3_lon))
      
         s = sprintf('%s %s %s %s %s %s', pinfo.var1_var, pinfo.var2_var, pinfo.var3_var, ...
                                              model, pinfo.fname, pinfo.ens_mem);
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
         outm{nlines+1} = sprintf('%s %s %s %s %s %s', pinfo.var1_var, ...
               pinfo.var2_var, pinfo.var3_var, model, pinfo.fname, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);
      
         set(objh(1),'interpreter','none')
      end
      legend boxoff

   otherwise

      error(sprintf('model %s not implemented yet', model))

end


%======================================================================                    
% Subfunctions                                                                             
%======================================================================                    
                                                                                           
function var = Get1Copy(fname, copyindex, var, lvlind, latind, lonind)                                            
% Gets a time-series of a single specified copy of a prognostic variable                   
% at a particular 3D location (level, lat, lon)                                            
if strcmp(var,'ps')                                                                  
   corner = [-1 copyindex        latind lonind];                 
   endpnt = [-1 copyindex        latind lonind];                 
else                                                                                       
   corner = [-1 copyindex lvlind latind lonind];                 
   endpnt = [-1 copyindex lvlind latind lonind];                 
end                                                                                        
var = getnc(fname, var, corner, endpnt);
