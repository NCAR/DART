% Check_pop_to_dart
%
% data.cal holds the starting time information
%

popdir   = '/ptmp/thoar/POP/job_40705';
popfile  = 'pop.r.x1A.00000102';
dartfile = 'assim_model_state_ud';

fname = sprintf('%s/%s',popdir,popfile);

S   = nc_varget(fname, 'SALT_CUR');
T   = nc_varget(fname, 'TEMP_CUR');
U   = nc_varget(fname, 'UVEL_CUR');
V   = nc_varget(fname, 'VVEL_CUR');
SSH = nc_varget(fname,'PSURF_CUR');

modelsize = prod(size(S)) + prod(size(T)) + ...
            prod(size(U)) + prod(size(V)) + prod(size(SSH));

[nx ny nz] = size(S);

iyear   = nc_attget(fname,nc_global,'iyear');
imonth  = nc_attget(fname,nc_global,'imonth');
iday    = nc_attget(fname,nc_global,'iday');
ihour   = nc_attget(fname,nc_global,'ihour');
iminute = nc_attget(fname,nc_global,'iminute');
isecond = nc_attget(fname,nc_global,'isecond');

fprintf('year  month  day  hour  minute  second %d %d %d %d %d %d\n',  ...
        iyear,imonth,iday,ihour,iminute,isecond);

% Get the dart equivalent

fname   = sprintf('%s/%s',popdir,dartfile);
fid     = fopen(fname,'rb','ieee-le');
trec1   = fread(fid,1,'int32');
seconds = fread(fid,1,'int32');
days    = fread(fid,1,'int32');
trecN   = fread(fid,1,'int32');

if (trec1 ~= trecN) 
   error('first record mismatch')
end

rec1    = fread(fid,        1,'int32');
datmat  = fread(fid,modelsize,'float64');
recN    = fread(fid,        1,'int32');
fclose(fid);

if (rec1 ~= recN) 
   error('second record mismatch')
end

disp(sprintf('time tag is days/seconds : %d %d',days,seconds))

offset = prod(size(S)); ind1 = 1; ind2 = offset;
myS    = reshape(datmat(ind1:ind2),size(S));
d      = myS - S;
disp(sprintf('S diffs are %f %f',min(d(:)),max(d(:))))

offset = prod(size(T)); ind1 = ind2+1; ind2 = ind1 + offset - 1;
myT    = reshape(datmat(ind1:ind2),size(T)); 
d      = myT - T;
disp(sprintf('T diffs are %f %f',min(d(:)),max(d(:))))

offset = prod(size(U)); ind1 = ind2+1; ind2 = ind1 + offset - 1;
myU    = reshape(datmat(ind1:ind2),size(U)); 
d      = myU - U;
disp(sprintf('U diffs are %f %f',min(d(:)),max(d(:))))

offset = prod(size(V)); ind1 = ind2+1; ind2 = ind1 + offset - 1;
myV    = reshape(datmat(ind1:ind2),size(V)); 
d      = myV - V;
disp(sprintf('V diffs are %f %f',min(d(:)),max(d(:))))

offset = prod(size(SSH)); ind1 = ind2+1; ind2 = ind1 + offset - 1;
mySSH  = reshape(datmat(ind1:ind2),size(SSH)); 
d      = mySSH - SSH;
disp(sprintf('SSH diffs are %f %f',min(d(:)),max(d(:))))

clear fid fname iyear imonth iday ihour iminute isecond
clear trec1 trecN rec1 recN offset ind1 ind2 days seconds
