function var_list()
%VAR_LIST:  variable list for the 9 var model

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

fprintf('-----------------------------------------\n');
var_type = strvcat('x (divergence)', 'y (height)', 'z (pressure)');

for i = 1:3
   for j = 1:3
      index = (i - 1) * 3 + j;
      fprintf('Variable %d is %10s %d\n', index, var_type(i, :), j);
   end
end

fprintf('-----------------------------------------\n');

