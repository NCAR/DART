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


