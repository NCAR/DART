function temp = compute_temperature( pres, theta, Cp, Rd, p0 )
%% FUNCTION compute_temperature - Computes temperature from potential temperature.
%
% Inputs:  
%	pres     = pressure, at mass pts
%       theta    = potential temperature, at mass pts 
%	Cp,Rd,p0 = c_p, dry gas constants, reference pressure
% Output:
%	temp     = temperature, at mass pts

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

kappa = Rd / Cp ;

temp = theta .* (pres ./ p0).^kappa ;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
