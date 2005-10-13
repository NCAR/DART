function temp = compute_temperature( pres, theta, Cp, Rd, p0 )
%
% Computes temperature from potential temperature.
%
% Inputs:  
%	pres     = pressure, at mass pts
%       theta    = potential temperature, at mass pts 
%	Cp,Rd,p0 = c_p, dry gas constants, reference pressure
% Output:
%	temp     = temperature, at mass pts

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.htm

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

kappa = Rd / Cp ;

temp = theta .* (pres ./ p0).^kappa ;
