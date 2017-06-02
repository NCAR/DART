function cov_factor = comp_cov_factor(z_in, c)
%% comp_cov_factor Gaspari Cohn cutoff, z_in is the distance while c is the cutoff

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

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

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
