function x = ChecknetCDFuse(fname)
% ChecknetCDFuse  checks to see if the netcdf pieces are available. 
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

fid = fopen(fname,'wt');
fprintf(fid,'%s \n',datestr(now));
fprintf(fid,'%s \n',pwd);
if (exist('getnc') ~= 2)
   fprintf(fid,'%s \n','No matlab netcdf operators, not using matlab.');
   fprintf(fid,'%d \n',-1);
   x = -1;
else
   fprintf(fid,'%s \n','Found matlab netcdf operators, using matlab.');
   fprintf(fid,'%d \n',0);
   x = 0;
end
fclose(fid);
