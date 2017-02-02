function height = compute_height( phi, g )
%% compute_height
%
% Inputs:
%	phi = (full) geopotential, at w pts
%	g   = gravitational acceleration
% Output:
%	height = height, at mass pts
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

height = ( phi(2:end,:,:) + phi(1:end-1,:,:) ) ./ (2*g) ;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
