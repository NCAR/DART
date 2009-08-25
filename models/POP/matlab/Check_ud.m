function [dart pop] = Check_ud(dartfile,popfile)
% Check_ud : check the conversion of a POP restart to a DART state vector.
% dartfile = 'assim_model_state_ud';
%  popfile = 'pop.r.nc' 
% x        = Check_ud(dartfile,popfile);
%
%
% dartfile = '/fs/image/home/thoar/SVN/DART/models/POP/work/assim_model_state_ud';
%  popfile = '/fs/image/home/thoar/SVN/DART/models/POP/work/cx3.dart.001.pop.r.0002-01-01-00000.nc';
% [dart pop] = Check_ud(dartfile,popfile);
%

% Read the original POP file values.
% The nc_varget() function returns the variables with the fastest 
% varying dimension on the right. This is opposite to the Fortran
% convention of the fastest varying dimension on the left ... so 
% one of the variables must be permuted in order to be compared.

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

Sdiff     = pop.S     - dart.S;     [min(    Sdiff(:)) max(    Sdiff(:))]
Tdiff     = pop.T     - dart.T;     [min(    Tdiff(:)) max(    Tdiff(:))]
Udiff     = pop.U     - dart.U;     [min(    Udiff(:)) max(    Udiff(:))]
Vdiff     = pop.V     - dart.V;     [min(    Vdiff(:)) max(    Vdiff(:))]
PSURFdiff = pop.PSURF - dart.PSURF; [min(PSURFdiff(:)) max(PSURFdiff(:))]


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
