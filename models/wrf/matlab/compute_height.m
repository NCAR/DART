function height = compute_height( phi, g )
%
% Inputs:
%	phi = (full) geopotential, at w pts
%	g   = gravitational acceleration
% Output:
%	height = height, at mass pts
%

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

height = ( phi(2:end,:,:) + phi(1:end-1,:,:) ) ./ (2*g) ;
