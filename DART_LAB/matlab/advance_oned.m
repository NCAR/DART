function x_new = advance_oned(x, alpha)

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

x_new = x + comp_dt(x, alpha);
end

%---------------------------------------------------

% Internal function comp_dt
function dx = comp_dt(x, alpha)

% Compute the time tendency; alpha controls nonlinearity
dx = x + alpha .* x .* abs(x);

end
