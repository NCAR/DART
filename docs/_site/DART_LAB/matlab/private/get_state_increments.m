function [state_incs, Rxy] = get_state_increments(state_ens, obs_ens, obs_incs)
%% get_state_increments Computes state increments given observation increments and
% the state and obs prior ensembles

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Compute state variance and covariance
covar = cov(state_ens, obs_ens);

Rxy = covar(1, 2);

state_incs = obs_incs * covar(1, 2) / covar(2, 2);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
