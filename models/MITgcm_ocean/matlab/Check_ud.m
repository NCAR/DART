function chunk = Check_ud(fname,dsize)
%% Check_ud a routine to check the (forcibly big-endian binary) files on a little-endian machine
%
% fname = 'assim_model_state_ud.0001';
% dsize = [256 225 40];
% x = Check_ud(fname,dsize);
%
% fname = '/ptmp/thoar/MITgcm/Run-2/advance_temp0/assim_model_state_ud';
% dsize = [256 225 40];
% xi = Check_ud(fname,dsize);
%
% fname = '/fs/image/home/thoar/SVN/DART/models/MITgcm_ocean/ibrahim/assim_model_state_ud';
% dsize = [256 225 40];
% xg = Check_ud(fname,dsize);
%
% Sdiff   = xi.S - xg.S; [min(Sdiff(:)) max(Sdiff(:))]
% Tdiff   = xi.T - xg.T; [min(Tdiff(:)) max(Tdiff(:))]
% Udiff   = xi.U - xg.U; [min(Udiff(:)) max(Udiff(:))]
% Vdiff   = xi.V - xg.V; [min(Vdiff(:)) max(Vdiff(:))]
% Etadiff = xi.Eta - xg.Eta; [min(Etadiff(:)) max(Etadiff(:))]

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

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

chunk.fname   = fname;
chunk.seconds = seconds;
chunk.days    = days;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
