function varid = SetCopyID(fname);
%% SetCopyID queries for the copy index for a set of ensemble members of a specific netCDF file.

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

metadata    = nc_varget(fname,'CopyMetaData');       % get all the metadata
copyindices = strmatch('ensemble member',metadata);  % find all 'ensemble member's

if ( isempty(copyindices) )
   fprintf('%s has no valid ensemble members\n',fname)
   disp('To be a valid ensemble member, the CopyMetaData for the member')
   disp('must start with the character string ''ensemble member''')
   disp('None of them in do in your file.')
   fprintf('%s claims to have %d copies\n',fname, num_copies)
   error('netcdf file has no ensemble members.')
end

ncopies    = length(copyindices);
def_copies = round([1*ncopies/4 , 2*ncopies/4 , 3*ncopies/4 ]);
def_string = sprintf(' %d ',def_copies);

disp('Enter any individual ensemble members IDs to plot.')
fprintf('2 4 13 (between 1 and %d)    ... or ...\n',ncopies)
disp('13                           ... or ... ')
disp('-1                           for none.')
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

