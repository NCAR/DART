%% ratio_test

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Plot the ratio of the new linear variance
% update factor for filter to the old linear in sd case.
% The factor gamma and the ratio of updated to prior
% regressed standard deviation, a, can go from 0 to 1.
% Look at ratio of new factor to old factor on this
% domain.

x(1:100) = 0.0;
y(1:100) = 0.0;
for i = 1:100
   x(i) = i / 100;
   y(i) = i / 100;
end

ratio(1:100, 1:100) = 1.0;
for i = 1 : 99,
   a = i / 100;
   for j = 1 : 100, 
      gamma = j / 100;
      gamma_p = (sqrt(gamma * a * a - gamma + 1) - 1) / (a - 1);
      ratio(i, j) = gamma_p / gamma;
   end
end

contourf(x, y, ratio, 10);
colorbar;
xlabel('regressed square root ratio, a');
ylabel('distance factor, gamma');
title('ratio of linear sd to linear variance updates');

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
