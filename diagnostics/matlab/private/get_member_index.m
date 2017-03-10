function member_index = get_member_index(fname, memstring)
%% GET_MEMBER_INDEX  Gets an index corresponding to MemberMetadata string
% Retrieves index associated with a given metadata string in
% file fname. If string is not found, a -1 is returned.
%
% Example:
% fname = 'filter_output.nc';
% memstring = 'ensemble member 5';
% member_index = get_member_index(fname, memstring);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

% Matlab seems to always need to transpose character variables.
metadata = ncread(fname,'MemberMetadata')';
[ens_size, ~]  = nc_dim_info(fname,'member');
[metalen, ~]   = nc_dim_info(fname,'metadatalength');

if( size(metadata,1) ~= ens_size || size(metadata,2) ~= metalen)
    error('%s from %s does not have the shape expected',memstring,fname)
end

nowhitecs = dewhite(memstring);

% Figure out which copy is the matching one
member_index = -1;
for i = 1:ens_size,

   % for matching -- we want to ignore whitespace -- find it & remove it
   nowhitemd = dewhite(metadata(i,:));

   if strcmp(nowhitemd , nowhitecs) == 1
      member_index = i;
   end
end

% Provide modest error support

if (member_index < 0)
   fprintf('ERROR: %s is not a valid metadata string for file %s\n', ...
                strtrim(memstring), fname)
   disp('valid metadata strings are: ')
   for i = 1:ens_size,
      fprintf('%s\n',deblank(metadata(i,:)))
   end
   error('Thats all. Start over')
end


function str2 = dewhite(str1)
%  function to remove ALL whitespace from a character string

str2 = str1(~isspace(str1));


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
