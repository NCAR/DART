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
rank = max(find(s_ens < squeeze(x))) + 1;;
if(isempty(rank))
   rank = 1;
end

end
