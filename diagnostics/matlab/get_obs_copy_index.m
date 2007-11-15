function indx = get_obs_copy_index(fid,copystring)
% get_obs_copy_index(fid,copystring)
% fid is the result of the 'netcdf' command 
% copystring is the desired copy 
%
% EXAMPLE:
% fid = netcdf('obs_seq.final.nc','nowrite');
% indx = get_obs_copy_index(fid,'bias');
%

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


copynames = fid{'copy'};
ncopies   = size(copynames,1);

indx      = -1; % Assume Failure

for i = 1:ncopies

   x = lower(deblank(copystring));
   y = lower(deblank(copynames(i,:)));

   if ( strcmp(x,y) )
      indx = i;
      break
   end 

end

if (indx < 1)
   error(sprintf('Did not find copy index for %s ',copystring))
end
