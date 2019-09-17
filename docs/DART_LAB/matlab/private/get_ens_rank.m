function [rank] = get_ens_rank(ens, x)
%% get_ens_rank

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

s_ens = sort(ens);
rank = find(s_ens < squeeze(x),1,'last') + 1;
if(isempty(rank))
   rank = 1;
end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
