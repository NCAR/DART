function temp = compute_temperature( pres, theta, Cp, Rd, p0 )
% FUNCTION compute_temperature - Computes temperature from potential temperature.
%
% Inputs:  
%	pres     = pressure, at mass pts
%       theta    = potential temperature, at mass pts 
%	Cp,Rd,p0 = c_p, dry gas constants, reference pressure
% Output:
%	temp     = temperature, at mass pts

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

kappa = Rd / Cp ;

temp = theta .* (pres ./ p0).^kappa ;
