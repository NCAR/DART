function copy_index = get_qc_index(fname, copystring, context)
%% get_qc_index  Gets an index corresponding to QCMetaData string
% Retrieves index associated with a given QCMetaData string in fname.
% If string does not match any, a fatal error is thrown. 
%
% Example:
% fname = 'obs_epoch_001.nc';
% copystring = 'Data QC';
% copy_index = get_qc_index(fname, copystring);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Need to get the QCMetaData (strings with the names)
% We then search the metadata for each copy for the appropriate copies.

errorstring = sprintf('\nERROR: "%s" is not a valid QCMetaData value for file %s\n', ...
              strtrim(copystring), fname);

if (nargin == 3)
   msgstring = sprintf('valid values for "%s" are', context);
else
   msgstring = 'valid values for QCMetaData are';
end

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

qc_meta_data    = ncread(fname,'QCMetaData')';
[num_copies, ~] = nc_dim_info(fname,'qc_copy');

% For a single copy, the size is nx1, for more k copies, it's kxn
if size(qc_meta_data, 2) == 1
   qc_meta_data = transpose(qc_meta_data);
end

nowhitecs = dewhite(copystring);

% Figure out which copy is the matching one
copy_index = -1;
for i = 1:num_copies,

   % for matching -- we want to ignore whitespace -- find it & remove it
   nowhitemd = dewhite(qc_meta_data(i,:));

   if strcmp(nowhitemd , nowhitecs) == 1
      fprintf('%s is QC copy %3i\n', copystring,i);
      copy_index = i;
   end
end

% Provide modest error support

if (copy_index < 0)
   for i = 1:num_copies,
      msgstring = sprintf('%s\n%s',msgstring,deblank(qc_meta_data(i,:)));
   end
   error(sprintf('%s\n%s',errorstring,msgstring))
end

function str2 = dewhite(str1)
% internal function to remove whitespace from a character string
i = find( ~isspace(str1) );
str2 = str1(i);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
