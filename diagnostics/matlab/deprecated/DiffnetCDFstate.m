function DiffnetCDFstate(file1,file2,outfile)
%% DiffnetCDFstate  checks to see if the netcdf pieces are available.
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
% file1   = 'test1/preassim.nc';
% file2   = 'test2/preassim.nc';
% outfile = 'test1_2.diff';
% DiffnetCDFstate(file1,file2,outfile)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%----------------------------------------------------------------------
% Check for the existence of the two input files.
%----------------------------------------------------------------------

if (exist(file1,'file') ~= 2), error('%s does not exist.',file1); end
if (exist(file2,'file') ~= 2), error('%s does not exist.',file2); end

%----------------------------------------------------------------------
% if the 'state' variable exists, read in the WHOLE THING.
%----------------------------------------------------------------------

if (    nc_isvar(file1,'state'))
   a = nc_varget(file1,'state');
else
   error('%s has no ''state'' variable.',file1)
end

if (    nc_isvar(file2,'state'))
   b = nc_varget(file2,'state');
else
   error('%s has no ''state'' variable.',file2)
end

% string them into 1D arrays and take the difference

c = b(:) - a(:);

%c(1) = c(1) + 0.1e-15;

a = min(c);
b = max(c);
clear c

% Write the min and max of the differences to a file

fid = fopen(outfile,'wt');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%s\n',pwd);
fprintf(fid,'%s\n',file1);
fprintf(fid,'%s\n',file2);

if ( (a == 0.0) && (b == 0.0) )
   fprintf(fid,'0 %e %e\n',a,b);
else
   fprintf(fid,'1 %e %e\n',a,b);
end

fclose(fid);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
