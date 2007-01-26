function DiffnetCDFstate(file1,file2,outfile)
% DiffnetCDFstate  checks to see if the netcdf pieces are available. 
% 
% DiffnetCDFstate(nc1, nc2, outfile);
%
% fname ... a filename that contains diagnostic information.
%           The last line contains three things:
%           1) an error code (0 == files identical)
%           2) the min of the difference
%           3) the max of the difference
%
% Example:
%
% file1   = 'test1/Prior_Diag.nc';
% file2   = 'test2/Prior_Diag.nc';
% outfile = 'test1_2.diff';
% DiffnetCDFstate(file1,file2,outfile)


% Data Assimilation Research Testbed -- DART
% Copyright 2004-2006, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% $Id: DiffnetCDFstate.m,v 1.1 2007/01/26 21:16:43 thoar Exp $

%----------------------------------------------------------------------
% Check for the existence of the two input files.
%----------------------------------------------------------------------

if (exist(file1,'file') ~= 2)
   error(sprintf('%s does not exist.',file1))
end
if (exist(file2,'file') ~= 2)
   error(sprintf('%s does not exist.',file2))
end

%----------------------------------------------------------------------
% if the 'state' variable exists, read in the WHOLE THING.
%----------------------------------------------------------------------

f1 = netcdf(file1);
f2 = netcdf(file2);

if (   VarExist(f1,'state'))
   a = getnc(file1,'state');
else
   error(sprintf('%s has no ''state'' variable.',file1))
end

if (   VarExist(f2,'state'))
   b = getnc(file2,'state');
else
   error(sprintf('%s has no ''state'' variable.',file2))
end

% string them into 1D arrays and take the difference

c = b(:) - a(:);

%c(1) = c(1) + 0.1e-15;

a = min(c);
b = max(c);
clear c f1 f2

% Write the min and max of the differences to a file

fid = fopen(outfile,'wt');
fprintf(fid,'%s \n',datestr(now));
fprintf(fid,'%s \n',pwd);
fprintf(fid,'%s \n',file1);
fprintf(fid,'%s \n',file2);

if ( (a == 0.0) && (b == 0.0) )
   fprintf(fid,'0 %e %e\n',a,b);
else
   fprintf(fid,'1 %e %e\n',a,b);
end

fclose(fid);


function x = VarExist(ncid,varname)

x = 0;   % false ... assumed not to exist.
variables = var(ncid);
for i=1:length(variables)
   if ( strmatch(name(variables{i}), varname) == 1 )
      x = 1;   % true ... variables exists in the netcdf file.
   end
end

