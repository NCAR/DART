function [dart pop] = Check_pop_to_dart(popfile,dartfile)
%% Check_pop_to_dart : check pop_to_dart.f90 ... the conversion of a POP restart to a DART state vector file.
%
%  popfile = 'pop.r.nc';
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

% The nc_varget() function returns the variables with the fastest 
% varying dimension on the right. This is opposite to the Fortran
% convention of the fastest varying dimension on the left ... so 
% one of the variables must be permuted in order to be compared.

S     = nc_varget(popfile,  'SALT_CUR'); pop.S     = permute(S,   [3 2 1]);
T     = nc_varget(popfile,  'TEMP_CUR'); pop.T     = permute(T,   [3 2 1]);
U     = nc_varget(popfile,  'UVEL_CUR'); pop.U     = permute(U,   [3 2 1]);
V     = nc_varget(popfile,  'VVEL_CUR'); pop.V     = permute(V,   [3 2 1]);
PSURF = nc_varget(popfile, 'PSURF_CUR'); pop.PSURF = permute(PSURF, [2 1]);

disp(sprintf('pop.PSURF min/max are %0.8g %0.8g',min(pop.PSURF(:)),max(pop.PSURF(:))))

[nx ny nz] = size(pop.U);
fprintf('vert dimension size is %d\n',nz)
fprintf('N-S  dimension size is %d\n',ny)
fprintf('E-W  dimension size is %d\n',nx)

modelsize = nx*ny*nz;

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
dart.S   = get_data(fid, [nx ny nz], 'float64');
dart.T   = get_data(fid, [nx ny nz], 'float64');
dart.U   = get_data(fid, [nx ny nz], 'float64');
dart.V   = get_data(fid, [nx ny nz], 'float64');
dart.SSH = get_data(fid, [nx ny   ], 'float64');
recN     = fread(fid,     1,  'int32');
fclose(fid);

fprintf(' shape of DART variables is %d \n',size(dart.S))

% The POP restart file has PSURF ... DART drags around SSH
% SSH = psurf/980.6;

dart.PSURF = dart.SSH * 980.6;

disp(sprintf('PSURF min/max are %0.8g %0.8g',min(dart.PSURF(:)),max(dart.PSURF(:))))
disp(sprintf('SSH   min/max are %0.8g %0.8g',min(dart.SSH(:)),max(dart.SSH(:))))

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
% Add something known to each state variable to check dart_to_pop.f90

S1     = dart.S     + 1.0;
T1     = dart.T     + 2.0;
U1     = dart.U     + 3.0;
V1     = dart.V     + 4.0;
SSH1   = dart.SSH + 5.0;

datvec = [S1(:); T1(:); U1(:); V1(:); SSH1(:)];
clear S1 T1 U1 V1 SSH1

fid     = fopen('test.ic','wb','ieee-le');
% Write the 'advance_to_time' FIRST
fwrite(fid,  trec1,'int32');
fwrite(fid,seconds,'int32');
fwrite(fid,   days+1,'int32');
fwrite(fid,  trecN,'int32');

% Write the 'model_state_time' (close to the data)
fwrite(fid,  trec1,'int32');
fwrite(fid,seconds,'int32');
fwrite(fid,   days,'int32');
fwrite(fid,  trecN,'int32');

% Write the (modified) model state ... 
fwrite(fid,   rec1, 'int32');
fwrite(fid, datvec, 'float64');
fwrite(fid,   recN, 'int32');
fclose(fid);

% That's all folks ...

function B = get_data(fid, shape, typestr)
A = fread(fid, prod(shape), typestr);
B = reshape(A, shape);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
