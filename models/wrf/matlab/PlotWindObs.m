%

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

for periods = 2:37;

   fname = sprintf('wind_vectors.%03d.dat',periods);
   ncname = 'obs_diag_output.nc';
   platform = 'SAT';
   level = -1;

   obs = plot_wind_vectors(fname,ncname,platform,level);

   disp('Pausing, hit any key to continue ...')
   pause

end

