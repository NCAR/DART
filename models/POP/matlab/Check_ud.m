function [dart pop] = Check_ud(popfile,dartfile)
%% Check_ud : check pop_to_dart.f90 ... the conversion of a POP restart to a DART state vector file.
%
%  popfile = 'pop.r.nc' 
% dartfile = 'dart_ics';
% x        = Check_pop_to_dart(popfile, dartfile);
%
%  popfile = '~DART/models/POP/work/cx3.dart.001.pop.r.0002-01-01-00000.nc';
% dartfile = '~DART/models/POP/work/perfect_ics';
% [dart pop] = Check_pop_to_dart(popfile, dartfile);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Read the original POP file values.
if (exist(popfile,'file') ~= 2)
   error('POP file %s does not exist.',popfile)
end
if (exist(dartfile,'file') ~= 2)
   error('DART file %s does not exist.',dartfile)
end

iyear   = nc_attget(popfile,nc_global,'iyear');
imonth  = nc_attget(popfile,nc_global,'imonth');
iday    = nc_attget(popfile,nc_global,'iday');
ihour   = nc_attget(popfile,nc_global,'ihour');
iminute = nc_attget(popfile,nc_global,'iminute');
isecond = nc_attget(popfile,nc_global,'isecond');

fprintf('POP year  month  day  hour  minute  second %d %d %d %d %d %d\n',  ...
        iyear,imonth,iday,ihour,iminute,isecond);

pop.S     = nc_varget(popfile,  'SALT_CUR');
pop.T     = nc_varget(popfile,  'TEMP_CUR');
pop.U     = nc_varget(popfile,  'UVEL_CUR');
pop.V     = nc_varget(popfile,  'VVEL_CUR');
pop.PSURF = nc_varget(popfile, 'PSURF_CUR');

[nz ny nx] = size(pop.U);
fprintf('vert dimension size is %d\n',nz)
fprintf('N-S  dimension size is %d\n',ny)
fprintf('E-W  dimension size is %d\n',nx)

modelsize = nz*ny*nx;

% filesize = S,T,U,V * (nx*ny*nz) + SSH * (nx*ny)
storage  = 8;
size2d   = nx*ny;
size3d   = nx*ny*nz;
n2ditems = 1*size2d;
n3ditems = 4*size3d;
rec1size = 4+(4+4)+4;  % time stamps ... 
rec2size = 4+(n3ditems*storage + n2ditems*storage)+4;

fsize    = rec1size + rec2size;
disp(sprintf('with a modelsize of %d the file size should be %d bytes', ...
     modelsize,fsize))

% Open and read timetag for state
fid     = fopen(dartfile,'rb','ieee-le');
trec1   = fread(fid,1,'int32');
seconds = fread(fid,1,'int32');
days    = fread(fid,1,'int32');
trecN   = fread(fid,1,'int32');

fprintf('need to know POP calendar for better comparison.\n', days,seconds);
fprintf('DART days seconds %d %d\n', days,seconds);

if (trec1 ~= trecN) 
   error('first record mismatch')
end

% Successively read state vector variables.
rec1     = fread(fid,     1,  'int32');
dart.S   = get_3D_permuted(fid, [nx ny nz], 'float64');
dart.T   = get_3D_permuted(fid, [nx ny nz], 'float64');
dart.U   = get_3D_permuted(fid, [nx ny nz], 'float64');
dart.V   = get_3D_permuted(fid, [nx ny nz], 'float64');
dart.SSH = get_2D_permuted(fid, [nx ny   ], 'float64');
recN     = fread(fid,     1,  'int32');
fclose(fid);

fprintf(' shape of DART variables is %d \n',size(dart.S))

% The POP restart file has PSURF ... DART drags around SSH
% SSH = psurf/980.6;

dart.PSURF = dart.SSH * 980.6;

if (rec1 ~= recN) 
   error('second record mismatch')
end

dart.dartfile = dartfile;
dart.seconds  = seconds;
dart.days     = days;

% Find the range of the mismatch

d = pop.S     - dart.S;     disp(sprintf('S     diffs are %0.8g %0.8g',min(d(:)),max(d(:))))
d = pop.T     - dart.T;     disp(sprintf('T     diffs are %0.8g %0.8g',min(d(:)),max(d(:))))
d = pop.U     - dart.U;     disp(sprintf('U     diffs are %0.8g %0.8g',min(d(:)),max(d(:))))
d = pop.V     - dart.V;     disp(sprintf('V     diffs are %0.8g %0.8g',min(d(:)),max(d(:))))
d = pop.PSURF - dart.PSURF; disp(sprintf('PSURF diffs are %0.8g %0.8g',min(d(:)),max(d(:))))

% As an added bonus, we create an 'assim_model_state_ic' file with an 
% advance-to-time one day in the future.

% Open and read timetag for state
fid     = fopen(dartfile,'rb','ieee-le');
trec1   = fread(fid,1,'int32');
seconds = fread(fid,1,'int32');
days    = fread(fid,1,'int32');
trecN   = fread(fid,1,'int32');

% read state vector variables.
rec1     = fread(fid,     1,  'int32');
datvec   = fread(fid, n3ditems+n2ditems, 'float64');
recN     = fread(fid,     1,  'int32');
fclose(fid);

% Open and write advance_to_time
fid     = fopen('test.ic','wb','ieee-le');
fwrite(fid,  trec1,'int32');
fwrite(fid,seconds,'int32');
fwrite(fid,   days,'int32');
fwrite(fid,  trecN,'int32');

fwrite(fid,  trec1,'int32');
fwrite(fid,seconds,'int32');
fwrite(fid,   days+1,'int32');
fwrite(fid,  trecN,'int32');

% read state vector variables.
fwrite(fid,   rec1, 'int32');
fwrite(fid, datvec, 'float64');
fwrite(fid,   recN, 'int32');
fclose(fid);



% The nc_varget() function returns the variables with the fastest 
% varying dimension on the right. This is opposite to the Fortran
% convention of the fastest varying dimension on the left ... so 
% one of the variables must be permuted in order to be compared.

function C = get_3D_permuted(fid, shape, typestr)
datasize = prod(shape);
A = fread(fid, prod(shape), typestr);
B = reshape(A, shape);
C = permute(B, [3 2 1]);

function C = get_2D_permuted(fid, shape, typestr)
datasize = prod(shape);
A = fread(fid, prod(shape), typestr);
B = reshape(A, shape);
C = permute(B, [2 1]);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
