function CompareObsDiag(file1,file2)
%% CompareObsDiag checks to see if two obs_diag netcdf files are identical.
%
% CompareObsDiag will subtract each numeric variable and print a 
% summary of the min/max differences for each variable.
% If there are a different number of 'copies', for example, you can
% specify which copies to compare in the call to get_hyperslab().
% At present, this is a manual edit within this file.
% This file is largely similar to Compare_netCDF_files.m but explicitly
% checks only the FIRST 21 copies and uses a snctools get_var instead
% of the native Matlab netcdf.getVar().
%
% Example:
%
% file1 = 'obs_diag_output.nc';
% file2 = 'obs_diag_output_proposed.nc';
% CompareObsDiag(file1,file2)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fid = fopen('CompareObsDiag_results.txt','a+');
fprintf(fid,'\n---\n');
fprintf(fid,'%s\n',datestr(now));

fprintf(fid,'Comparing %s\n',file1);
fprintf(fid,'with      %s\n',file2);
CompareTwoFiles(file1,file2,fid)

fclose(fid);

%%---------------------------------------------------------------------
%  Display min and max of the difference of two numeric variables 
%  from two separate files.
%----------------------------------------------------------------------

function CompareTwoFiles(file1, file2, fid)

if (exist(file1,'file') ~= 2), error('%s does not exist.',file1); end
if (exist(file2,'file') ~= 2), error('%s does not exist.',file2); end

f1info     = ncinfo(file1);
nvariables = length(f1info.Variables);
ncid1      = netcdf.open(file1,'NOWRITE');
ncid2      = netcdf.open(file2,'NOWRITE');

% Loop over all the variables in the first file ... the variables
% can be in any order in either file. I do not check to see if
% there are variables present in the second file that are not in
% the first file.

for ivar = 1:nvariables

   varid1 = ivar-1; % Because netCDF tools use C-like indices.

   [varname1, xtype1, ~, numatts1] = netcdf.inqVar(  ncid1,varid1);
   varid2                          = netcdf.inqVarID(ncid2,varname1);
   [varname2, xtype2, ~, numatts2] = netcdf.inqVar(  ncid2,varid2);

   FillValue1 = has_att(ncid1, varid1, numatts1, '_FillValue');
   FillValue2 = has_att(ncid2, varid2, numatts2, '_FillValue');

   switch (f1info.Variables(ivar).Datatype)
      case {'char','character'}
         fprintf(fid,'%32s is a %s variable.\n',varname1,f1info.Variables(ivar).Datatype);
         fprintf(    '%32s is a %s variable.\n',varname1,f1info.Variables(ivar).Datatype);
      otherwise
         
         bob   = get_hyperslab('fname',file1,'varname',varname1,'copy1',1,'copycount',21);
         data1 = bob(:);
         bob   = get_hyperslab('fname',file2,'varname',varname1,'copy1',1,'copycount',21);
         data2 = bob(:);
        
         bob     = data1 - data2;
         datamin = min(bob);
         datamax = max(bob);

         fprintf(fid,'%32s has min/max differences of %g %g\n',varname1,datamin,datamax);
         if ((datamin ~= 0.0) && (datamax ~= 0.0))
             fprintf('%32s has min/max differences of %g %g\n',varname1,datamin,datamax);
         end
   end

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

function data = my_getVar(ncid,varid,FillValue,start,count)

data = netcdf.getVar(ncid, varid, start, count);
if ( ~ isempty(FillValue) )
   data(data == FillValue) = NaN;
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
