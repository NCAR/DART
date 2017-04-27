function copy_index = get_qc_index(fname, copystring)
%% get_qc_index  Gets an index corresponding to QC meta_data string
% Retrieves index associated with a given QC meta_data string in
% file fname. If string is not found in meta_data list, a -1 is returned.
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

% Need to get the QC MetaData (strings with the names)
% We then search the metadata for each copy for the appropriate copies.

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
   fprintf('ERROR: %s is not a valid metadata string for file %s\n', ...
                copystring, fname)
   disp('valid metadata strings are: ')
   for i = 1:num_copies,
      disp(sprintf('%s',deblank(qc_meta_data(i,:))))
   end
   error('please try again')
end

function str2 = dewhite(str1)
% internal function to remove whitespace from a character string
i = find( ~isspace(str1) );
str2 = str1(i);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
