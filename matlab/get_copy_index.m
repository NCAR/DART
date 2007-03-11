function copy_index = get_copy_index(fname, copystring)
%GET_COPY_INDEX  Gets an index corresponding to copy meta_data string
% Retrieves index associated with a given copy meta_data string in 
% file fname. If string is not found in meta_data list, a -1 is returned.
%
% Example:
% fname = 'Prior_Diag.nc';
% copystring = 'ensemble member 5';
% copy_index = get_copy_index(fname, copystring);

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

f = netcdf(fname);
num_copies = ncsize(f{'copy'}); % determine # of ensemble members
close(f)

% Need to get the Copy MetaData (strings with the names)
% We then search the metadata for each copy for the appropriate copies.
copy_meta_data = getnc(fname, 'CopyMetaData');

% For a single copy, the size is nx1, for more k copies, it's kxn
if size(copy_meta_data, 2) == 1
   copy_meta_data = transpose(copy_meta_data);
end

nowhitecs = dewhite(copystring);

% Figure out which copy is the matching one
copy_index = -1;
for i = 1:num_copies,

   % for matching -- we want to ignore whitespace -- find it & remove it
   nowhitemd = dewhite(copy_meta_data(i,:));

   if strcmp(nowhitemd , nowhitecs) == 1
      fprintf('%s is copy %3i\n', copystring,i);
      copy_index = i;
   end
end

% Provide modest error support

if (copy_index < 0)
   disp(sprintf('ERROR: %s is not a valid metadata string for file %s', ...
                copystring, fname))
   disp('valid metadata strings are: ')
   for i = 1:num_copies,
      disp(sprintf('%s',deblank(copy_meta_data(i,:))))
   end
   error('please try again')
end

function str2 = dewhite(str1)
% internal function to remove whitespace from a character string
i = find( ~isspace(str1) );
str2 = str1(i);
