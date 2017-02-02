function x = ChecknetCDFuse(fname)
%% ChecknetCDFuse  checks to see if the netcdf pieces are available.
%
% x = ChecknetCDFuse(fname);
%
% fname ... a filename that contains diagnostic information.
%           The last line contains a '0' (everything OK),
%           or a nonzero character ... matlab netcdf is not available.
%
% x  ... return code ... 0 == normal termination (OK)
%                     nonzero == bad, matlab netcdf not available.
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

fid = fopen(fname,'wt');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%s\n',pwd);
if (exist('nc_varget') ~= 2)
   fprintf(fid,'%s \n','No matlab snctools, not using matlab.');
   fprintf(fid,'%d \n',-1);
   x = -1;
else
   fprintf(fid,'%s \n','Found matlab snctools, using matlab.');
   fprintf(fid,'%d \n',0);
   x = 0;
end
fclose(fid);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
