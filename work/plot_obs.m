%

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next four lines automatically updated by CVS, do not edit>
% $Source$ 
% $Revision$ 
% $Date$ 
% $Author$ 

load obs_def;
x = obs_def(:, 1);
y = obs_def(:, 2);
plot(x, y, '*')
axis([0 360 -90 90])
