function PlotVarVarCorrel( pinfo )
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.
% The correlation is done across ensemble members.
%
% PlotVarVarCorrel is intended to be called by 'plot_var_var_correl'
%
% USAGE: PlotVarVarCorrel(fname, base_var_index, base_tme, state_var_index)
%
% fname
% base_var_index    index of one of the state variables to correlate  
% base_tme         index of the time of interest (timestep)
% state_var_index   index of the other state variable of interest.
%
% Example  (lorenz 63 model with 1000 timesteps)
%%--------------------------------------------------------
% pinfo.fname = 'Posterior_Diag.nc';
% pinfo.base_var_index    = 2;
% pinfo.base_tme         = 500;
% pinfo.state_var_index   = 1;
% PlotVarVarCorrel( pinfo )

% TJH Wed Jul  2 09:52:18 MDT 2003

if (exist(pinfo.fname) ~= 2), error(sprintf('%s does not exist.',pinfo.fname)), end

% disp(sprintf('PlotVarVarCorrel: fname is %s',pinfo.fname))
% disp(sprintf('PlotVarVarCorrel: base_var_index is %d',pinfo.base_var_index))
% disp(sprintf('PlotVarVarCorrel: base_tme      is %d',pinfo.base_tme))

% Get some file-specific information.
f = netcdf(pinfo.fname,'nowrite');
model      = f.model(:);
timeunits  = f{'time'}.units(:);
num_vars   = ncsize(f('StateVariable')); % determine # of state variables
num_times  = ncsize(f('time')); % determine # of output times
num_copies = ncsize(f('copy')); % determine # of ensemble members
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
          model, pinfo.base_var, pinfo.base_tme, pinfo.base_lvl, ...
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
      plot([pinfo.base_tme pinfo.base_tme],[ -1 1 ],'k:', ...
           [ax(1)         ax(2)],[  0 0 ],'k:')


   otherwise

      disp(sprintf('PlotVarVarCorrel: num_vars is %d',num_vars))
      
      % The Base Variable Index must be a valid state variable
      
      if ( pinfo.base_var_index > num_vars )
         disp( sprintf('%s only has %d state variables', pinfo.fname, num_vars))
         error(sprintf('you wanted variable # %d ', pinfo.base_var_index))
      end
      
      % The Time must be within range also.
      
      if ( pinfo.base_tme > num_times )
         disp( sprintf('%s only has %d output times', pinfo.fname, num_times))
         error(sprintf('you wanted time # %d ', pinfo.base_tme))
      end
      
      % The State Variable Index must also be a valid state variable
      
      if ( pinfo.base_var_index > num_vars )
         disp( sprintf('%s only has %d state variables', pinfo.fname, num_vars))
         error(sprintf('you wanted variable # %d ', pinfo.base_var_index))
      end
      
      % Get 'standard' ensemble series 
       base_var = get_ens_series(pinfo.fname,  pinfo.base_var_index);
      state_var = get_ens_series(pinfo.fname, pinfo.state_var_index);
      
      % perform a single correlation
      correl = ens_correl(base_var, pinfo.base_tme, state_var);
      
      clf; plot(correl);
      
      s1 = sprintf('%s Correlation of state variable %d, T = %d, with variable %d', ...
               model, pinfo.base_var_index, pinfo.base_tme, pinfo.state_var_index);
      s2 = sprintf('%d ensemble members -- %s', num_copies-2,pinfo.fname); 
      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel('time (timestep #)')
      ylabel('correlation')
      
      % call out the time index in question, and put a corr==0 reference line.
      ax = axis;
      hold on;
      plot([pinfo.base_tme pinfo.base_tme],[ -1 1 ],'k:', ...
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

