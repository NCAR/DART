function PlotVarVarCorrel( pinfo )
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.
% The correlation is done across ensemble members.
%
% PlotVarVarCorrel is intended to be called by 'plot_var_var_correl'
%
% USAGE: PlotVarVarCorrel( pinfo )
%
% pinfo is a structure with the following necessary components:
% fname
% base_var_index    index of one of the state variables to correlate  
% base_time         index of the time of interest (timestep)
% state_var_index   index of the other state variable of interest.
%
% Example  (lorenz 63 model with 1000 timesteps)
%%--------------------------------------------------------
% pinfo.fname = 'Posterior_Diag.nc';
% pinfo.base_var          = 'state';
% pinfo.base_var_index    = 2;
% pinfo.base_time         = 500;
% pinfo.state_var         = 'state';
% pinfo.state_var_index   = 1;
% PlotVarVarCorrel( pinfo )

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

if (exist(pinfo.fname) ~= 2), error(sprintf('%s does not exist.',pinfo.fname)), end

% Get some file-specific information.
f = netcdf(pinfo.fname,'nowrite');
model      = f.model(:);
timeunits  = f{'time'}.units(:);
var_atts   = dim(f{pinfo.base_var}); % cell array of dimensions for the var
num_times  = length(var_atts{1});    % determine # of output times
num_copies = length(var_atts{2});    % determine # of ensemble members
num_vars   = length(var_atts{3});    % determine # of state variables (of this type)
close(f);

switch lower(model)

   case 'fms_bgrid'

      clf;

      base_mem = GetEns( pinfo.fname, pinfo.base_var, ...
                    pinfo.base_lvlind, pinfo.base_latind, pinfo.base_lonind );
      comp_mem = GetEns( pinfo.fname, pinfo.comp_var, ...
                    pinfo.comp_lvlind, pinfo.comp_latind, pinfo.comp_lonind );

      correl = ens_correl(base_mem, pinfo.base_tmeind, comp_mem);
      times  = getnc(pinfo.fname,'time');

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
         plot(times,correl);

      s1 = sprintf('%s Correlation of ''%s'', T = %d, lvl = %d, lat = %.2f, lon=%.2f', ...
          model, pinfo.base_var, pinfo.base_time, pinfo.base_lvl, ...
          pinfo.base_lat, pinfo.base_lon);

      s2 = sprintf('with ''%s'', lvl = %d, lat = %.2f, lon= %.2f, %d ensemble members -- %s', ...
          pinfo.comp_var, pinfo.comp_lvl, pinfo.comp_lat, pinfo.comp_lon, ...
          num_copies-2,pinfo.fname); 

      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel(sprintf('time (%s) %d timesteps',timeunits, num_times))
      ylabel('correlation')
      
      % call out the time index in question, and put a corr==0 reference line.
      ax = axis;
      hold on;
      plot([pinfo.base_time pinfo.base_time],[ -1 1 ],'k:', ...
           [ax(1)         ax(2)],[  0 0 ],'k:')


   otherwise

      disp(sprintf('PlotVarVarCorrel: num_vars is %d',num_vars))
      
      % The Base Variable Index must be a valid state variable
      
      if ( pinfo.base_var_index > num_vars )
         disp( sprintf('%s only has %d state variables', pinfo.fname, num_vars))
         error(sprintf('you wanted variable # %d ', pinfo.base_var_index))
      end
      
      % The Time must be within range also.
      
      if ( pinfo.base_time > num_times )
         disp( sprintf('%s only has %d output times', pinfo.fname, num_times))
         error(sprintf('you wanted time # %d ', pinfo.base_time))
      end
      
      % The State Variable Index must also be a valid state variable
      
      if ( pinfo.base_var_index > num_vars )
         disp( sprintf('%s only has %d state variables', pinfo.fname, num_vars))
         error(sprintf('you wanted variable # %d ', pinfo.base_var_index))
      end
      
      % Get 'standard' ensemble series 
       base_var = get_ens_series(pinfo.fname, pinfo.base_var,  pinfo.base_var_index);
      state_var = get_ens_series(pinfo.fname, pinfo.state_var, pinfo.state_var_index);
      
      % perform a single correlation
      correl = ens_correl(base_var, pinfo.base_time, state_var);
      
      clf; plot(correl);
      
      s1 = sprintf('%s Correlation of variable %s %d, T = %d, with variable %s %d', ...
               model, pinfo.base_var, pinfo.base_var_index, pinfo.base_time, ...
                      pinfo.state_var, pinfo.state_var_index);
      s2 = sprintf('%d ensemble members -- %s', num_copies-2,pinfo.fname); 
      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel('time (timestep #)')
      ylabel('correlation')
      
      % call out the time index in question, and put a corr==0 reference line.
      ax = axis;
      hold on;
      plot([pinfo.base_time pinfo.base_time],[ -1 1 ],'k:', ...
           [ax(1)         ax(2)],[  0 0 ],'k:')
      
      %axis(ax)
end


%======================================================================
% Subfunctions
%======================================================================

function var = GetEns( fname, var, lvlind, latind, lonind)
% Gets a time-series of all copies of a prognostic variable 
% at a particular 3D location (level, lat, lon).
% Determining just the ensemble members (and not mean, spread ...)
% is the hard part.

% find which are actual ensemble members
metadata    = getnc(fname,'CopyMetaData');           % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's

if ( isempty(copyindices) )
   disp(sprintf('%s has no valid ensemble members',fname))
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   disp(sprintf('%s claims to have %d copies',fname, num_copies))
   error('netcdf file has no ensemble members.')
end
ens_num     = length(copyindices);

% Get all ensemble members, just return desired ones.
if strcmp(var,'ps')
   corner = [-1 -1        latind lonind];
   endpnt = [-1 -1        latind lonind];
else
   corner = [-1 -1 lvlind latind lonind];
   endpnt = [-1 -1 lvlind latind lonind];
end
bob = getnc(fname, var, corner, endpnt); % 'bob' is only 2D 
var = bob(:,copyindices);


function PlotLocator(pinfo)
   plot(pinfo.base_lon, pinfo.base_lat,'pg','MarkerSize',12,'MarkerFaceColor','g');
   axis([0 360 -90 90])
   worldmap
   axis image
   grid on

