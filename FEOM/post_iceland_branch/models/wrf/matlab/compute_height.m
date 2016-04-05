function height = compute_height( phi, g )
%
% Inputs:
%	phi = (full) geopotential, at w pts
%	g   = gravitational acceleration
% Output:
%	height = height, at mass pts
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.htm

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

height = ( phi(2:end,:,:) + phi(1:end-1,:,:) ) ./ (2*g) ;
