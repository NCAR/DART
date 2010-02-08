function PlotJeffCorrel( pinfo )
%% Plots exploratory correlation plots. Don't use without talking to J. Anderson.
%
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.
% The correlation is done across ensemble members.
%
% PlotVarVarCorrel is intended to be called by 'plot_var_var_correl'
%
% USAGE: PlotVarVarCorrel(fname, base_var_index, base_time, state_var_index)
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

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist(pinfo.fname) ~= 2), error('%s does not exist.',pinfo.fname), end

% Get some file-specific information.

model      = nc_attget(pinfo.fname, nc_global, 'model');
timeunits  = nc_attget(pinfo.fname, 'time',    'units');
num_times  = dim_length(pinfo.fname,'time');
num_copies = dim_length(pinfo.fname,'copy');

switch lower(model)

   case {'fms_bgrid','pe2lyr','mitgcm_ocean'}

      clf;

      base_mem = GetEns( pinfo.fname, pinfo.base_var, ...
                    pinfo.base_lvlind, pinfo.base_latind, pinfo.base_lonind );
      comp_mem = GetEns( pinfo.fname, pinfo.comp_var, ...
                    pinfo.comp_lvlind, pinfo.comp_latind, pinfo.comp_lonind );
      nmembers = size(comp_mem,2);

      correl = jeff_correl(base_mem, pinfo.base_tmeind, comp_mem);
      times  = nc_varget(pinfo.fname,'time');

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
         plot(times,correl);

      s1 = sprintf('%s Correlation of ''%s'', T = %d, lvl = %d, lat = %.2f, lon=%.2f', ...
          model, pinfo.base_var, pinfo.base_time, pinfo.base_lvl, ...
          pinfo.base_lat, pinfo.base_lon);

      s2 = sprintf('with ''%s'', lvl = %d, lat = %.2f, lon= %.2f, %d ensemble members -- %s', ...
          pinfo.comp_var, pinfo.comp_lvl, pinfo.comp_lat, pinfo.comp_lon, ...
          nmembers, pinfo.fname); 

      title({s1,s2},'interpreter','none','fontweight','bold')
      xlabel(sprintf('time (%s) %d timesteps',timeunits, num_times))
      ylabel('correlation')
      
      % call out the time index in question, and put a corr==0 reference line.
      ax = axis;
      hold on;
      plot([pinfo.base_time pinfo.base_time],[ -1 1 ],'k:', ...
           [ax(1)         ax(2)],[  0 0 ],'k:')


   otherwise

      num_vars   = dim_length(pinfo.fname,'StateVariable');
      fprintf('PlotJeffCorrel: num_vars is %d\n',num_vars)
      
      % The Base Variable Index must be a valid state variable
      
      if ( pinfo.base_var_index > num_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, num_vars)
         error('you wanted variable # %d ', pinfo.base_var_index)
      end
      
      % The Time must be within range also.
      
      if ( pinfo.base_time > num_times )
         fprintf('%s only has %d output times\n', pinfo.fname, num_times)
         error('you wanted time # %d ', pinfo.base_time)
      end
      
      % The State Variable Index must also be a valid state variable
      
      if ( pinfo.base_var_index > num_vars )
         fprintf('%s only has %d state variables\n', pinfo.fname, num_vars)
         error('you wanted variable # %d ', pinfo.base_var_index)
      end
      
      % Get 'standard' ensemble series 
       base_var = get_ens_series(pinfo.fname, pinfo.base_var,  pinfo.base_var_index);
      state_var = get_ens_series(pinfo.fname, pinfo.state_var, pinfo.state_var_index);
      nmembers  = size(state_var,2);

      % perform a single correlation
      correl = jeff_correl(base_var, pinfo.base_time, state_var);
      
      clf; plot(correl);
      
      s1 = sprintf('%s Correlation of variable %s %d, T = %d, with variable %s %d', ...
               model, pinfo.base_var, pinfo.base_var_index, pinfo.base_time, ...
                      pinfo.state_var, pinfo.state_var_index);
      s2 = sprintf('%d ensemble members -- %s', nmembers, pinfo.fname); 
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



function x = dim_length(fname,dimname)
bob = nc_getdiminfo(fname,dimname);
x   = bob.Length;



function var = GetEns( fname, var, lvlind, latind, lonind)
% Gets a time-series of all copies of a prognostic variable 
% at a particular 3D location (level, lat, lon).
% Determining just the ensemble members (and not mean, spread ...)
% is the hard part.

% find which are actual ensemble members
metadata    = nc_varget(fname,'CopyMetaData');       % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'member's

if ( isempty(copyindices) )
   fprintf('%s has no valid ensemble members\n',fname)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',fname, num_copies)
   error('netcdf file has no ensemble members.')
end
ens_num     = length(copyindices);

% Get all ensemble members, just return desired ones.
myinfo.diagn_file = fname;
myinfo.levelindex = lvlind;
myinfo.latindex   = latind;
myinfo.lonindex   = lonind;
[start, count]    = GetNCindices(myinfo,'diagn',var);

bob = nc_varget(fname, var, start, count); % 'bob' is only 2D 
var = bob(:,copyindices);



function PlotLocator(pinfo)
   plot(pinfo.base_lon, pinfo.base_lat,'pb','MarkerSize',12,'MarkerFaceColor','b');
   hold on;
   plot(pinfo.comp_lon, pinfo.comp_lat,'pr','MarkerSize',12,'MarkerFaceColor','r');
   hold off;
   axis([0 360 -90 90]);
   worldmap;
   axis image
   grid on

