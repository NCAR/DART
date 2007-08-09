function ostruct = CombineStructs(struct1,struct2);
% CombineStructs   all components of both stuctures are combined into one structure.
%
% EXAMPLE:
%
%

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ~( isstruct(struct1) & isstruct(struct2) )
   error('both arguments must be structures')
end

% We'll just copy one input structure to the output
% and append all unique fields from #2 into the output.
% If the field exists in both, and is different, we're in trouble.

ostruct = struct1;
fields  = fieldnames(struct2);

for i=1:length(fields)

   ostruct = setfield(ostruct, fields{i}, getfield(struct2,fields{i}));

end
