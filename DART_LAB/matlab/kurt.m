function k = kurt(vals)
% computes the kurtosis of the given input array
%
% based on the second formula on this web page:
% http://www.ats.ucla.edu/stat/mult_pkg/faq/general/kurtosis.htm  

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% compute mean  - WARNING: this assumes the array is laid out
% so the number of items is the second dimension, e.g. [1, 4]

nvals = size(vals, 2);
m = sum(vals) / nvals;

% compute s2 and s4 
del = vals - m;
s2 = sum(del .* del);                % element-by-element *
s4 = sum(del .* del .* del .* del);  % ditto
    
% compute m2 and m4 
m2 = s2 / nvals;
m4 = s4 / nvals;

% finally, the kurtosis value.  this is the version of kurtosis
% that does NOT subtract 3.0 from the result.
k = (m4 ./ (m2 .* m2));

end 
