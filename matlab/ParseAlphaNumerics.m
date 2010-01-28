function [vrbl, vrbl_inds] = ParseAlphaNumeric(IDstring)
% ParseAlphaNumerics -  extricates a variable name from subsequent IDs 
% str1 = ' X 1 3 4 89'
% [alpha, numerics] = ParseAlphaNumerics(str1)
% alpha = 'X'
% numerics = [1 3 4 89];

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

inds       = find(IDstring == ',');     % find all commas
IDstring(inds) = ' ';
words      = strread(IDstring,'%s');
nwords     = length(words);
vrbl       = words{1};

for i = 2:nwords
   vrbl_inds(i-1) = sscanf(words{i},'%d');
end
