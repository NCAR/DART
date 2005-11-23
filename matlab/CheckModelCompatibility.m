function CheckModelCompatibility(truth_file,diagn_file);
% CheckModelCompatibility   tries to ensure that two netcdf files can be compared.
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

if ( exist(truth_file) ~= 2 )
   error(sprintf('(truth_file) %s does not exist.',truth_file))
end
if ( exist(diagn_file) ~= 2 )
   error(sprintf('(diagn_file) %s does not exist.',diagn_file))
end

% Get some information from the truth_file
% dimensions are ft(xxxx), variables are ft{xxxx}
ft = netcdf(truth_file);
tmodel      = ft.model(:);
tnum_vars   = ncsize(ft('StateVariable')); % determine # of state variables
tnum_copies = ncsize(ft('copy')); % determine # of ensemble members
tnum_times  = ncsize(ft('time')); % determine # of output times
close(ft); 

if (isempty(tmodel)) 
   error(sprintf('%s has no ''model'' global attribute.',truth_file))
end
if (prod(size(tnum_vars)) > 1 ) 
   error(sprintf('%s has no ''StateVariable'' dimension.',truth_file))
end
if (prod(size(tnum_copies)) > 1 ) 
   error(sprintf('%s has no ''copy'' dimension.',truth_file))
end
if (prod(size(tnum_times)) > 1 ) 
   error(sprintf('%s has no ''time'' dimension.',truth_file))
end

% Get some information from the diagn_file
fd = netcdf(diagn_file);
dmodel      = fd.model(:);
dnum_vars   = ncsize(fd('StateVariable')); % determine model size
dnum_copies = ncsize(fd('copy')); % determine # of ensemble members
dnum_times  = ncsize(fd('time')); % determine # of output times
close(fd); 

if (isempty(dmodel)) 
   error(sprintf('%s has no ''model'' global attribute.',diagn_file))
end
if (prod(size(dnum_vars)) > 1 ) 
   error(sprintf('%s has no ''StateVariable'' dimension.',diagn_file))
end
if (prod(size(dnum_copies)) > 1 ) 
   error(sprintf('%s has no ''copy'' dimension.',diagn_file))
end
if (prod(size(dnum_times)) > 1 ) 
   error(sprintf('%s has no ''time'' dimension.',diagn_file))
end

% rudimentary bulletproofing
if (strcmp(tmodel,dmodel) ~= 1)
   disp(sprintf('%s has model %s ',truth_file,tmodel))
   disp(sprintf('%s has model %s ',diagn_file,dmodel))
   error('no No NO ... models must be the same')
end
if (tnum_vars ~= dnum_vars)
   disp(sprintf('%s has %d state variables',truth_file,tnum_vars))
   disp(sprintf('%s has %d state variables',diagn_file,dnum_vars))
   error('no No NO ... both files must have same number of state variables.')
end
if (tnum_times ~= dnum_times)
   disp(sprintf('%s has %d timesteps',truth_file,tnum_times))
   disp(sprintf('%s has %d timesteps',diagn_file,dnum_times))
   error('ugh ... both files must have same number of timesteps.')
end
