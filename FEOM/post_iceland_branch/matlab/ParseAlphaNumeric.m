function [vrbl, vrbl_inds] = ParseAlphaNumeric(IDstring)
% ParseAlphaNumeric    local function 
% to extricate a variable name from subsequent IDs 
% str1 = ' X 1 3 4 89'
% [alpha, numerics] = ParseAlphaNumeric(str1)
% alpha = 'X'
% numerics = [1 3 4 89];

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

str1       = upper(IDstring);       % convert to uppercase
inds       = find(str1 ~= ' ');     % find all non-blanks
vrbl       = str1(inds(1));         % use first non-blank char

inds       = find(str1 == vrbl);
str1(inds) = ' ';                   % remove variable from string
vrbl_inds  = sscanf(str1,'%d');
vrbl_inds  = reshape(vrbl_inds,[1,length(vrbl_inds)]);
