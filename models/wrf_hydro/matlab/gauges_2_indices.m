function obs_ind_id = gauges_2_indices(gauges)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id: gauges_2_indices.m $

n_links = length(gauges);

k = 0;
for i = 1:n_links
    ob_id = gauges(i, :);
    
    if sum( isspace(ob_id) ) ~= 10
        k = k + 1;
        obs_ind_id(k, :) = [i, str2double(ob_id)]; %#ok
    end
end

% <next few lines under version control, do not edit>
% $URL: $
% $Revision: $
% $Date: $
