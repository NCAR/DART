function PlotPhaseSpace( pinfo )
% PlotPhaseSpace: Plots trajectories of 3 variables for any ensemble member
%
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
% var1name     name of state variable to plot as 'X'
% var1ind      ID   of state variable to plot as 'X'
% var2name     name of state variable to plot as 'Y'
% var2ind      ID   of state variable to plot as 'Y'
% var3name     name of state variable to plot as 'Z'
% var3ind      ID   of state variable to plot as 'Z'
% ltype        line type (see 'help plot' for details)
%
% Example 1   ( 9 variable model )
%%--------------------------------------------------------
% pinfo.fname      = 'True_State.nc';
% pinfo.ens_mem    = 'true state';    % true state only has 1 ens mem ...
% pinfo.var1name   = 'state';         % 9var netCDF has only 1 flavor variable 
% pinfo.var2name   = 'state';    
% pinfo.var3name   = 'state';    
% pinfo.var1ind    = 3;    
% pinfo.var2ind    = 6;
% pinfo.var3ind    = 7;
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

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( exist(pinfo.fname) ~= 2 ), error(sprintf('file %s does not exist.',pinfo.fname)), end

% Get some information for the 'X' variable
f            = netcdf(pinfo.fname);
model        = f.model(:);
X.var_atts   = dim(f{pinfo.var1name});  % cell array of dimensions for the var
X.num_times  = length(X.var_atts{1});   % determine # of output times
X.num_copies = length(X.var_atts{2});   % # of ensemble members
X.num_vars   = length(X.var_atts{3});   % dimension of desired variable

% Get some information for the 'Y' variable
Y.var_atts   = dim(f{pinfo.var2name});  % cell array of dimensions for the var
Y.num_times  = length(Y.var_atts{1});   % determine # of output times
Y.num_copies = length(Y.var_atts{2});   % # of ensemble members
Y.num_vars   = length(Y.var_atts{3});   % dimension of desired variable

if ( isfield( pinfo,'var3name') )
   % Get some information for the 'Z' variable
   Z.var_atts   = dim(f{pinfo.var3name});  % cell array of dimensions for the var
   Z.num_times  = length(Z.var_atts{1});   % determine # of output times
   Z.num_copies = length(Z.var_atts{2});   % # of ensemble members
   Z.num_vars   = length(Z.var_atts{3});   % dimension of desired variable
   close(f);
end

switch lower(model)

   case {'ikeda'}   % only two state variables

      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);  % errors out if no ens_mem 
      
      x = get_var_series(pinfo.fname, pinfo.var1name, ens_mem_id, pinfo.var1ind);
      y = get_var_series(pinfo.fname, pinfo.var2name, ens_mem_id, pinfo.var2ind);

      h = plot(x,y,pinfo.ltype);
      axis image

      % datmin = min(min(x,y));
      % datmax = max(max(x,y));
      % axis([datmin datmax datmin datmax])
      
      % If there is no legend, we are assuming this is the first plot on this
      % axis. We need to add a title, legend, axes labels, etc.
      [legh, objh, outh, outm] = legend;
      if (isempty(legh))
      
         % title(sprintf('%s ensemble member %d of %s',model,ens_mem_id,pinfo.fname), ...
         %    'interpreter','none','fontweight','bold')
         title('The Ikeda phase space','fontweight','bold')
         xlabel(sprintf('%s variable # %d',pinfo.var1name, pinfo.var1ind))
         ylabel(sprintf('%s variable # %d',pinfo.var2name, pinfo.var2ind))
      
         s = sprintf('%d %d %s %s %s', pinfo.var1ind, pinfo.var2ind, ...
                     pinfo.fname, model, pinfo.ens_mem);
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
         outm{nlines+1} = sprintf('%d %d %s %s %s', ...
              pinfo.var1ind, pinfo.var2ind, pinfo.fname, ...
              model, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);
 
         set(objh(1:nlines+1),'interpreter','none')
      end
      legend boxoff

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','simple_advection'} 

      BulletProof(pinfo,X,Y,Z)          % rudimentary bulletproofing

      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);  % errors out if no ens_mem 
      
      x = get_var_series(pinfo.fname, pinfo.var1name, ens_mem_id, pinfo.var1ind);
      y = get_var_series(pinfo.fname, pinfo.var2name, ens_mem_id, pinfo.var2ind);
      z = get_var_series(pinfo.fname, pinfo.var3name, ens_mem_id, pinfo.var3ind);

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
         xlabel(sprintf('%s variable # %d',pinfo.var1name, pinfo.var1ind))
         ylabel(sprintf('%s variable # %d',pinfo.var2name, pinfo.var2ind))
         zlabel(sprintf('%s variable # %d',pinfo.var3name, pinfo.var3ind))
      
         s = sprintf('%d %d %d %s %s %s', pinfo.var1ind, pinfo.var2ind, ...
                     pinfo.var3ind, pinfo.fname, model, pinfo.ens_mem);
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
              pinfo.var1ind, pinfo.var2ind, pinfo.var3ind, pinfo.fname, ...
              model, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);
 
         set(objh(1:nlines+1),'interpreter','none')
      end
      legend boxoff

   case 'fms_bgrid'

      disp(sprintf('PlotPhaseSpace'))
      pinfo

      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);   % errors out if no ens_mem 
      
      x = Get1Copy(pinfo.fname, ens_mem_id, pinfo.var1name, ...
                  pinfo.var1_lvlind, pinfo.var1_latind, pinfo.var1_lonind);
      y = Get1Copy(pinfo.fname, ens_mem_id, pinfo.var2name, ...
                  pinfo.var2_lvlind, pinfo.var2_latind, pinfo.var2_lonind);
      z = Get1Copy(pinfo.fname, ens_mem_id, pinfo.var3name, ...
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
               pinfo.var1name, pinfo.var1_lvl, pinfo.var1_lat, pinfo.var1_lon))
         ylabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var2name, pinfo.var2_lvl, pinfo.var2_lat, pinfo.var2_lon))
         zlabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var3name, pinfo.var3_lvl, pinfo.var3_lat, pinfo.var3_lon))
      
         s = sprintf('%s %s %s %s %s %s', pinfo.var1name, pinfo.var2name, pinfo.var3name, ...
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
         outm{nlines+1} = sprintf('%s %s %s %s %s %s', pinfo.var1name, ...
               pinfo.var2name, pinfo.var3name, model, pinfo.fname, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);
      
         set(objh(1:nlines+1),'interpreter','none')
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



function BulletProof(pinfo,X,Y,Z)

if ( (pinfo.var1ind > X.num_vars) | (pinfo.var1ind < 1) ) 
   disp(sprintf('\n%s has %d %s variables',pinfo.fname,X.num_vars,pinfo.var1name))
   disp(sprintf('%d  <= ''var1'' <= %d',1,X.num_vars))
   error(sprintf('var1 (%d) out of range',pinfo.var1ind))
end

if ( (pinfo.var2ind > Y.num_vars) | (pinfo.var2ind < 1) ) 
   disp(sprintf('\n%s has %d %s variables',pinfo.fname,Y.num_vars,pinfo.var2name))
   disp(sprintf('%d  <= ''var2'' <= %d',1,Y.num_vars))
   error(sprintf('var2 (%d) out of range',pinfo.var2ind))
end

if ( (pinfo.var3ind > Z.num_vars) | (pinfo.var3ind < 1) ) 
   disp(sprintf('\n%s has %d %s variables',pinfo.fname,Z.num_vars,pinfo.var3name))
   disp(sprintf('%d  <= ''var3'' <= %d',1,Z.num_vars))
   error(sprintf('var3 (%d) out of range',pinfo.var3ind))
end
