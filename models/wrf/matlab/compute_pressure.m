function pres = compute_pressure( mu, dnw, phi, theta, qv, Rd,Rv,gamma,p0 )
%
% CS version of Ryan's m-file.  Computes pressure from gas law as in WRF.
%
% Inputs:  
%	mu    =   (full) dry hydrostatic surf. press (2d)
%	dnw   =   intervals between w levels, at mass pts
%	phi   =   (full) geopotential, at w pts
%       theta =   (full) potential temperature, at mass pts 
%       qv    =   water-vapor mixing ratio, at mass pts 
%	Rd,Rv,gamma,p0
%	      =   dry and moist gas constants, c_p/c_v, reference pressure
% Output:
%	pres  =   (full) pressure, at mass pts
%
% See wrf subroutine calc_p_rho_phi.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.htm

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

 % alpha = inverse density of dry air, computed from hydrostatic relation
 % for the dry air:  - dphi/deta = -mu alpha.

 phi_eta =  ...
  ( phi(2:end,:,:) - phi(1:end-1,:,:) ) ./ repmat( dnw(:), [1 size(mu)] );
 alpha = -phi_eta ./ repmat( reshape( mu, [1 size(mu)] ), [ length(dnw) 1 1 ] );

 % Gas law: 
 pres = p0 * ( Rd * theta .* (1 + (Rd/Rv)*qv ) ./ ( p0 * alpha ) ).^gamma ;
