function rho = compute_density( mu, dnw, phi )
% FUNCTION rho : Computes density of dry air, computed from the
%    hydrostatic relation for dry air: - dphi/deta = -mu alpha.
%
% Inputs:  
%	mu  = (full) dry hydrostatic surf. press (2d)
%	dnw = intervals between w levels, at mass pts
%	phi = (full) geopotential, at w pts
% Output:
%	rho  =  density, at mass pts
%
% See wrf subroutine calc_p_rho_phi.

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

phi_eta = ( phi(2:end,:,:) - phi(1:end-1,:,:) ) ./ repmat( dnw(:), [1 size(mu)] );

rho = -repmat( reshape( mu, [1 size(mu)] ), [ length(dnw) 1 1 ] ) ./ phi_eta ;
