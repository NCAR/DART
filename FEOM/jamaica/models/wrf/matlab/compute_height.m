function height = compute_height( phi, g )
%
% Inputs:
%	phi = (full) geopotential, at w pts
%	g   = gravitational acceleration
% Output:
%	height = height, at mass pts
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

height = ( phi(2:end,:,:) + phi(1:end-1,:,:) ) ./ (2*g) ;
