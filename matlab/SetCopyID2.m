function varid = SetCopyID2(fname);
% SetCopyID queries for the copy index for a set of ensemble members of a specific netCDF file.

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

metadata = nc_varget(fname,'CopyMetaData');           % get all the metadata
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
   else                              % we want these
      varid = ensmems;
   end
end

