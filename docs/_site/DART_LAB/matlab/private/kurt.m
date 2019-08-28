function k = kurt(vals)
% computes the kurtosis of the given input array
%
% based on the second formula on this web page:
% http://www.ats.ucla.edu/stat/mult_pkg/faq/general/kurtosis.htm

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% array of diffs from mean
del = vals - mean(vals);

% compute the square and 4th power of the diffs from mean
m2 = mean(del .^ 2);
m4 = mean(del .^ 4);

% compute the kurtosis value.  this is the version
% of the kurtosis formula that is not nonbiased and
% does not subtract 3.0 from the result.
k = (m4 ./ (m2 .^ 2));

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
