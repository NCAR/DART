function varid = SetCopyID(fname);
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

metadata    = getnc(fname,'CopyMetaData');           % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'ensemble member's

if ( isempty(copyindices) )
   disp(sprintf('%s has no valid ensemble members',fname))
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   disp(sprintf('%s claims to have %d copies',fname, num_copies))
   error('netcdf file has no ensemble members.')
end


ncopies    = length(copyindices);
def_copies = round([1*ncopies/4 , 2*ncopies/4 , 3*ncopies/4 ]);
def_string = sprintf(' %d ',def_copies);

disp('Enter any individual ensemble members IDs to plot.')
disp(sprintf('2 4 13 (between 1 and %d)    ... or ... ',ncopies))
disp(sprintf('13                           ... or ... ',ncopies))
disp(sprintf('-1                           for none.',ncopies))
disp('(no intervening syntax required)');
IDstring = input(sprintf('<cr> for %s\n',def_string),'s');

if isempty(IDstring)                 % take the default
   ensmems = def_copies;
   varid = zeros(1,length(ensmems));
   for i = 1:length(ensmems),
      copystring = sprintf('ensemble member %d',ensmems(i));
      varid(i) = get_copy_index(fname,copystring);
   end
else 
   ensmems = sscanf(IDstring,'%d');  % convert text to numbers
   if ( ensmems(1) < 0 )             % dont want any
      varid = [];
   else                              % we want these
      varid = zeros(1,length(ensmems));
      for i = 1:length(ensmems),
         copystring = sprintf('ensemble member %d',ensmems(i));
         varid(i) = get_copy_index(fname,copystring);
      end
   end
end
