function [state_incs] = get_state_increments(state_ens, obs_ens, obs_incs)
% Computes state increments given observation increments and
% the state and obs prior ensembles

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Compute state variance and covariance
covar = cov(state_ens, obs_ens);


state_incs = obs_incs * covar(1, 2) / covar(2, 2);

