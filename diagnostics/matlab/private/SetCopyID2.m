function [varid copymetadata] = SetCopyID2(fname);
%% SetCopyID2 queries for the copy indices in a specific netCDF file.
%  pays no attention to what the copy metadata means.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

metadata = ncread(fname,'MemberMetadata')'; % get all the metadata
ncopies  = size(metadata,1);

if ( ncopies < 2 )
   fprintf('%s has no valid ensemble members\n',fname)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',fname, num_copies)
   error('netcdf file has no ensemble members.')
elseif (ncopies == 2)
   def_copies = [1 2];
else
   def_copies = round([1*ncopies/3 , 2*ncopies/3]);
end

disp('Enter any (and all) IDs to plot.')
for i = 1:ncopies
   fprintf('ID %2d  is  %s\n',i,deblank(metadata(i,:)))
end
fprintf('-1                           for none.\n',ncopies)
disp('(no intervening syntax required)');

IDstring = input(sprintf('<cr> for %d %d \n',def_copies(1),def_copies(2)),'s');

if isempty(IDstring)                 % take the default
   varid = def_copies;
else
   ensmems = sscanf(IDstring,'%d');  % convert text to numbers
   if ( ensmems(1) < 0 )             % dont want any
      varid = [];
      copymetadata = [];
   else                              % we want these
      varid = ensmems;
   end
end

for i = 1:length(varid)
   copymetadata{i} = metadata(varid(i),:);
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
