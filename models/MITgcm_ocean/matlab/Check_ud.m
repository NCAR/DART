function chunk = Check_ud(fname,dsize)
%
% fname = 'assim_model_state_ud.0001';
% dsize = [256 225 40];
%
% x = Check_ud(fname,dsize);

modelsize = prod(dsize);

% S,T,U,V * dsize + SSH (dsize(1:2))

fsize = 4+(4+4)+4 + 4+(4*(modelsize*8)+(dsize(1)*dsize(2))*8)+4;
disp(sprintf('with a modelsize of %d the file size should be %d bytes', ...
     modelsize,fsize))

% Open and read timetag for state
fid     = fopen(fname,'rb','ieee-be');
trec1   = fread(fid,1,'int32');
seconds = fread(fid,1,'int32');
days    = fread(fid,1,'int32');
trecN   = fread(fid,1,'int32');

if (trec1 ~= trecN) 
   error('first record mismatch')
end

size3d = prod(dsize);
size2d = prod(dsize(1:2));
% Successively read state vector variables.
rec1     = fread(fid,        1,'int32');
datmat   = fread(fid,size3d,'float64'); chunk.S   = reshape(datmat,dsize);
datmat   = fread(fid,size3d,'float64'); chunk.T   = reshape(datmat,dsize);
datmat   = fread(fid,size3d,'float64'); chunk.U   = reshape(datmat,dsize);
datmat   = fread(fid,size3d,'float64'); chunk.V   = reshape(datmat,dsize);
datmat   = fread(fid,size2d,'float64'); chunk.Eta = reshape(datmat,dsize(1:2));
recN     = fread(fid,        1,'int32');
fclose(fid);

if (rec1 ~= recN) 
   error('second record mismatch')
end

chunk.fname = fname;
