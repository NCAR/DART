% plot_jeff_correl.m
% Plots time series of correlation between a given variable at a given
% time and another variable at all times in an ensemble time sequence.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% $Source$
% $Revision$
% $Date$

if (exist('fname') ~=1)
   disp('Input name of file;')
   fname = input('<cr> for Prior_Diag.nc\n','s');
   if isempty(fname)
      fname = 'Prior_Diag.nc';
   end
else
   % check to make sure they are using a file with some ensemble members.

   ft         = netcdf(fname);                                                     
   num_copies = ncsize(ft('copy')); % determine # of ensemble members              
   close(ft)

   if (num_copies <= 3) 
      error(sprintf('Sorry -- %s does not have enough ensemble members to correlate.',fname))
   end 
end 

vars = CheckModel(fname);   % also gets default values for this file.
pinfo.fname = fname;

switch lower(vars.model)

   case {'9var','lorenz_63','lorenz_96'}

      inputstring = input(sprintf('Input index for base variable (between %d and %d)  ', ...
                               vars.min_state_var, vars.max_state_var), 's');
      pinfo.base_var_index = str2num(deblank(inputstring));

      inputstring = input(sprintf('Input time index for base point (between 1 and %d)  ', ...
                               vars.time_series_length),'s');
      pinfo.base_tme = str2num(deblank(inputstring));

      inputstring = input(sprintf('Input variable index for correlation (between %d and %d)  ', ...
                               vars.min_state_var, vars.max_state_var), 's');
      pinfo.state_var_index = str2num(deblank(inputstring));

      disp(sprintf('Using diagnostic file %s',fname))
      disp(sprintf('Correlating state variable %d at time %d with state variable %d.', ...
              pinfo.base_var_index, pinfo.base_tme, pinfo.state_var_index))

   case 'fms_bgrid'

      pinfo = GetBgridInfo(fname, 'PlotVarVarCorrel');

      pinfo                            % just echo stuff for posterity.

   otherwise

      error(sprintf('model %s not implemented yet', vars.model))

end

% could/should check input for valid range, etc.

PlotJeffCorrel( pinfo )
clear vars inputstring
