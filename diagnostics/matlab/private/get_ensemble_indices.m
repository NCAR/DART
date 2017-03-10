function [ens_size, ens_indices] = get_ensemble_indices(fname)
%% DART:GET_ENSEMBLE_INDICES  returns the number of ensemble members in the file and their 'copy' indices. 
%
% Example:
% fname = 'filter_output.nc';
% [ens_size, ens_indices] = get_ensemble_indices(fname);
%
% Example to return just the size ...
% [ens_size, ~] = get_ensemble_indices(fname);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

[ens_size, ~] = nc_dim_info(fname,'member');

ens_indices = 1:ens_size;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
