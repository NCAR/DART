function [rank] = get_ens_rank(ens, x)

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

s_ens = sort(ens);

for i = 1 : size(ens, 2)
   if(x < s_ens(1, i))
      rank = i;
      return;
   end
end

rank = size(ens, 2) + 1;

end
