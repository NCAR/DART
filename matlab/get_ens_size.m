function ens_size = get_ens_size(fname)
%% GET_ENS_SIZE  returns the number of ensemble members in the netCDF file.
%
% Example:
% fname = 'Prior_Diag.nc';
% ens_size = get_ens_size(fname);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

copy_meta_data = ncread(fname,'MemberMetadata')';
atts           = nc_getdiminfo(fname,'copy');
num_copies     = atts.Length;

% For a single copy, the size is nx1, for more k copies, it's kxn
if size(copy_meta_data, 2) == 1
   copy_meta_data = transpose(copy_meta_data);
end

%% Figure out which copy is the matching one - ignoring whitespace
ens_size = 0;
for i = 1:num_copies,

   % for matching -- we want to ignore whitespace -- find it & remove it
   nowhitemd = dewhite(copy_meta_data(i,:));

   if strmatch('ensemblemember', nowhitemd) == 1
      ens_size = ens_size + 1;
   end
end

% Provide modest error support

if (ens_size == 0)
   fprintf('WARNING: unable to find any valid ensemble members in %s\n', fname)
   disp('valid metadata strings are: ')
   for i = 1:num_copies,
      fprintf('%s\n',deblank(copy_meta_data(i,:)))
   end
end

function str2 = dewhite(str1)
% internal function to remove whitespace from a character string
i = find( ~isspace(str1) );
str2 = str1(i);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
