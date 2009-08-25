function cov_factor = comp_cov_factor(z_in, c)
% Gaspari Cohn cutoff, z_in is the distance while c is the cutoff

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

z = abs(z_in);

if( z >= c*2.0)
   cov_factor = 0;
elseif( z <= c )
   r = z / c;
   cov_factor = ((( -0.25*r +0.5)*r +0.625)*r -5.0/3.0)*r^2 + 1.0;
else
   r = z / c;
   cov_factor = ((((r/12 -0.5)*r +0.625)*r +5.0/3.0)*r -5.0)*r + 4.0 - 2.0 / (3.0 * r);
end

