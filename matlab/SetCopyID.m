function varid = SetCopyID(fname)
%% SetCopyID queries for the copy index for a set of ensemble members of a specific netCDF file.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

[ncopies,copyindices] = get_ensemble_indices(fname);

if ( isempty(copyindices) )
   fprintf('%s has no ensemble members\n',fname)
   disp('To be a valid ensemble member, the MemberMetadata for the member')
   disp('must start with the character string ''ensemble member''')
   error('%s has no ensemble members.',fname)
end

def_copies  = round([1*ncopies/4 , 2*ncopies/4 , 3*ncopies/4 ]);
def_string  = sprintf(' %d ',def_copies);

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
      varid(i) = get_member_index(fname,copystring);
   end
else
   ensmems = sscanf(IDstring,'%d');  % convert text to numbers
   if ( ensmems(1) < 0 )             % dont want any
      varid = [];
   else                              % we want these
      varid = zeros(1,length(ensmems));
      for i = 1:length(ensmems),
         copystring = sprintf('ensemble member %d',ensmems(i));
         varid(i) = get_member_index(fname,copystring);
      end
   end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
