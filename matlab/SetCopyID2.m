function varid = SetCopyID2(fname);
% SetCopyID   Interactively determines the copy index for a set of 
%             ensemble members -- valid for the specific netCDF file.

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

metadata = getnc(fname,'CopyMetaData');           % get all the metadata
ncopies  = size(metadata,1);

if ( ncopies < 2 )
   disp(sprintf('%s has no valid ensemble members',fname))
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   disp(sprintf('%s claims to have %d copies',fname, num_copies))
   error('netcdf file has no ensemble members.')
elseif (ncopies == 2)
   def_copies = [1 2];
end

def_copies = round([1*ncopies/3 , 2*ncopies/3]);

disp('Enter any (and all) IDs to plot.')
for i = 1:ncopies
   disp(sprintf('ID %2d  is  %s',i,deblank(metadata(i,:))))
end
disp(sprintf('-1                           for none.',ncopies))
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
