%% DART:t_script

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

%x = -5:0.01:5;
%y = normcdf(x);

% Ensemble size is n
n = 20;
n1 = n + 1;


x = 1/n1 : 1/n1: n/n1;
y = norminv(x);

% Adjust the variance to get unit variance
y = (y - mean(y)) .* sqrt(1.0 / var(y)) + mean(y);

[mean(y) var(y) skewness(y) kurtosis(y)]


