function corr = jeff_correl(base_ens, comp_ens)
%% jeff_correl  Computes time-evolution of the correlation of a variable to another.
%
%  base_ens    is the ensemble state at many times for a given location
%  comp_ens   is another ensemble state at many times and a location.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

[base_nT, base_ens_size] = size(base_ens);
[comp_nT, comp_ens_size] = size(comp_ens);

if (base_ens_size ~= comp_ens_size)
   error('base ensemble size (%s) is not same as the comparison ensemble size (%s)',base_ens_size, comp_ens_size)
end
if (base_nT ~= comp_nT)
   error('base time series length (%s) is not same as the comparison time series length (%s)',base_nT, comp_nT)
end

corr    = zeros(base_nT,1);
corr(:) = NaN;

for i = 1:base_nT
   x = corrcoef(base_ens(i,:), comp_ens(i, :));
   corr(i) = x(1, 2);
end 


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
