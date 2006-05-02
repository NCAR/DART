function x = ChecknetCDFuse(fname)

fid = fopen(fname,'wt');
fprintf(fid,'%s \n',datestr(now));
fprintf(fid,'%s \n',pwd);
if (exist('getnc') ~= 2)
   fprintf(fid,'%s \n','No matlab netcdf operators, not using matlab.');
   fprintf(fid,'%d \n',-1);
else
   fprintf(fid,'%s \n','Found matlab netcdf operators, using matlab.');
   fprintf(fid,'%d \n',0);
end
fclose(fid);
