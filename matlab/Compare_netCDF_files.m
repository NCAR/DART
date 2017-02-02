function Compare_netCDF_files(file1,file2)
%% Compare_netCDF_files checks to see if two netcdf files are identical.
% An output text file of the summary of the min/max (differences) for each
% variable is produced.
%
% If called with a single argument, Compare_netCDF_files will print a
% summary of the min and max of every numeric variable in the file.
% If the file was the result of an 'ncdiff', for example, these should
% all be zero.
%
% If called with two arguments, Compare_netCDF_files will subtract each numeric
% variable and print a summary of the min/max differences for each variable to
% a text file. Only those variables with non-zero differences are output to 
% the command window.
%
% Example:
%
% file1   = '/glade/scratch/thoar/clm_2day1hop/run/clm_2day1hop.clm2_0004.r.2000-02-02-00000.nc';
% file2   = '/glade/scratch/thoar/clm_2day2hop/run/clm_2day2hop.clm2_0004.r.2000-02-02-00000.nc';
% Compare_netCDF_files(file1,file2)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fid = fopen('Compare_netCDF_files_results.txt','a+');
fprintf(fid,'\n---\n');
fprintf(fid,'%s\n',datestr(now));

if (nargin == 1)
   fprintf(fid,'Checking single file %s\n',file1);
   CheckSingleFile(file1,fid)
elseif(nargin == 2)
   fprintf(fid,'Comparing %s\n',file1);
   fprintf(fid,'with      %s\n',file2);
   CompareTwoFiles(file1,file2,fid)
else
   error('wrong number of arguments')
end

fclose(fid);

%%---------------------------------------------------------------------
%  Display min and max of all the numeric variables in a single file.
%----------------------------------------------------------------------

function CheckSingleFile(file1,fid)

if (exist(file1,'file') ~= 2), error('%s does not exist.',file1); end

f1info     = ncinfo(file1);
nvariables = length(f1info.Variables);
ncid1      = netcdf.open(file1,'NOWRITE');

for ivar = 1:nvariables

   varid = ivar-1; % Because netCDF tools use C-like indices.

   [varname, ~, ~, numatts] = netcdf.inqVar(ncid1,varid);
   FillValue = has_att(ncid1, varid, numatts, '_FillValue');

   switch (f1info.Variables(ivar).Datatype)
      case {'char','character'}
         fprintf(fid,'%32s is a %s variable.\n', varname, ...
                            f1info.Variables(ivar).Datatype);
         fprintf(    '%32s is a %s variable.\n', varname, ...
                            f1info.Variables(ivar).Datatype);
      otherwise
         data    = my_getVar(ncid1, varid, FillValue);
         bob     = data(:);
         datamin = min(bob);
         datamax = max(bob);
         fprintf(fid,'%32s has min/max of %g %g\n', varname, datamin, datamax);
         fprintf(    '%32s has min/max of %g %g\n', varname, datamin, datamax);
   end % switch
end

netcdf.close(ncid1)

%%---------------------------------------------------------------------
%  Display min and max of the difference of two numeric variables
%  from two separate files.
%----------------------------------------------------------------------

function CompareTwoFiles(file1, file2, fid)

if (exist(file1,'file') ~= 2), error('%s does not exist.',file1); end
if (exist(file2,'file') ~= 2), error('%s does not exist.',file2); end

f1info     = ncinfo(file1);
f2info     = ncinfo(file2);

if (length(f1info.Variables) ~= length(f2info.Variables) )
   fprintf('WARNING : the files have different numbers of variables.\n')
   fprintf('WARNING : only comparing the variables that exist in\n')
   fprintf('%s\n',file1)
end

nvariables = length(f1info.Variables);
ncid1      = netcdf.open(file1,'NOWRITE');
ncid2      = netcdf.open(file2,'NOWRITE');

% grab the variable name from one file and 
% search the second file for matching variable.

for ivar = 1:nvariables

   varid1 = ivar-1; % Because netCDF tools use C-like indices.

   [varname1, ~, ~, numatts1] = netcdf.inqVar(ncid1,varid1);
   varid2                     = netcdf.inqVarID(ncid2,varname1);
   [       ~, ~, ~, numatts2] = netcdf.inqVar(ncid2,varid2);
   
   FillValue1 = has_att(ncid1, varid1, numatts1, '_FillValue');
   FillValue2 = has_att(ncid2, varid2, numatts2, '_FillValue');

   switch (f1info.Variables(ivar).Datatype)
      case {'char','character'}
         fprintf(fid,'%32s is a %s variable.\n', varname1, ...
                            f1info.Variables(ivar).Datatype);
         fprintf(    '%32s is a %s variable.\n', varname1, ...
                            f1info.Variables(ivar).Datatype);
      otherwise
         bob     = my_getVar(ncid1, varid1, FillValue1);
         data1   = bob(:);
         bob     = my_getVar(ncid2, varid2, FillValue2);
         data2   = bob(:);

         if(numel(data1) ~= numel(data2))
            error('%32s has different number of elements %d and %d.\n', ...
                     varname1,numel(data1),numel(data2));
         end

         bob     = data1 - data2;
         datamin = min(bob);
         datamax = max(bob);

         fprintf(fid,'%32s has min/max differences of %g %g\n', ...
                     varname1,datamin,datamax);
         if ((datamin ~= 0.0) || (datamax ~= 0.0))
             fprintf('%32s has min/max differences of %g %g\n', ...
                     varname1,datamin,datamax);
         end
   end % switch
end

netcdf.close(ncid1)
netcdf.close(ncid2)

%%---------------------------------------------------------------------
%  helper functions
%----------------------------------------------------------------------

function data = has_att(ncid,varid,numatts,attstring)

data = [];

for iatt = 1:numatts
   attid = iatt - 1;
   attname = netcdf.inqAttName(ncid, varid, attid);

   switch( attname )
      case (attstring)
         data = netcdf.getAtt(ncid, varid, attstring);
      otherwise
   end
end

function data = my_getVar(ncid,varid,FillValue)

data = netcdf.getVar(ncid, varid);
if ( ~ isempty(FillValue) )
   data(data == FillValue) = NaN;
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
