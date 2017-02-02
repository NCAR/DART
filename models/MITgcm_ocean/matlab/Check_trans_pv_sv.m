function Check_trans_pv_sv
%% Check_trans_pv_sv()  routine to verify the packing/unpacking of the DART state vector.
% there are a lot of hardwired values ... 

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

dt = 900;
startDate_1=19960101;
startDate_2=60000;
tbase = datenum(1996,1,1,6,0,0);

timestep = 40992; % data70levels
timestep =  1344; % data40levels
timestep =    96; % real data  96*900 = 86400 = 1day

mitdir = '/ptmp/thoar/MITgcm/adv1day_assim/advance_temp0';

S   = rdmds(sprintf(  '%s/S.%10.10d',mitdir,timestep));
T   = rdmds(sprintf(  '%s/T.%10.10d',mitdir,timestep));
U   = rdmds(sprintf(  '%s/U.%10.10d',mitdir,timestep));
V   = rdmds(sprintf(  '%s/V.%10.10d',mitdir,timestep));
SSH = rdmds(sprintf('%s/Eta.%10.10d',mitdir,timestep));

modelsize = prod(size(S)) + prod(size(T)) + ...
            prod(size(U)) + prod(size(V)) + prod(size(SSH));

fsize = 4+4+4+4 + 4+modelsize*8+4;
disp(sprintf('with a modelsize of %d the file size should be %d bytes', ...
     modelsize,fsize))

toffset = (timestep*dt)/86400;

[nx ny nz] = size(S);

fname = 'assim_model_state_ud';
fname = '/ptmp/thoar/MITgcm/adv1day_assim/assim_model_state_ud.0020';
fid   = fopen(fname,'rb','ieee-be');

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

clear datmat myS myT myU myV mySSH d

%----------------------------------------------------------------------
% The perfect_ics file contains precisely two records.
% record 1 is two integers defining the valid time of the state vector.
% record 2 is the state vector.
% This is exactly the same as the 'assim_model_state_ud' file.
% (except they normally are for different times ...)
%----------------------------------------------------------------------

tbase = datenum(1601,1,1,0,0,0);  % this is zero in the DART(gregorian) world
t_one = datenum(1996,1,1,0,0,0);  % valid time of the 'gom' files.

toffset = t_one - tbase;
days    = floor(toffset);
seconds = round(toffset - days)*86400;

disp(sprintf('Creating a ''perfect_ics'' file with elements of shape %d %d %d',nx,ny,nz))

disp(sprintf('prod(size(S)) is %d',  prod(size(S))))
disp(sprintf('prod(size(T)) is %d',  prod(size(T))))
disp(sprintf('prod(size(U)) is %d',  prod(size(U))))
disp(sprintf('prod(size(V)) is %d',  prod(size(V))))
disp(sprintf('prod(size(SSH)) is %d',prod(size(SSH))))

nitems = prod(size(S))+ prod(size(T)) + prod(size(U)) + prod(size(V)) + prod(size(SSH));
disp(sprintf('total restart size should %d',nitems*8 + 8 + 16))

fid     = fopen('gom_S_199601.bin','rb','ieee-be');
[Sics,count] = fread(fid,prod(size(S)),'float32');
fclose(fid);
if (count ~= prod(size(S)))
   error(sprintf('S record length wrong %d %d',count,prod(size(S))))
end

fid     = fopen('gom_T_199601.bin','rb','ieee-be');
[Tics,count]    = fread(fid,prod(size(T)),'float32');
fclose(fid);
if (count ~= prod(size(T)))
   error(sprintf('T record length wrong %d %d',count,prod(size(T))))
end

fid     = fopen('gom_U_199601.bin','rb','ieee-be');
[Uics,count]    = fread(fid,prod(size(U)),'float32');
fclose(fid);
if (count ~= prod(size(U)))
   error(sprintf('U record length wrong %d %d',count,prod(size(U))))
end

fid     = fopen('gom_V_199601.bin','rb','ieee-be');
[Vics,count]    = fread(fid,prod(size(V)),'float32');
fclose(fid);
if (count ~= prod(size(V)))
   error(sprintf('V record length wrong %d %d',count,prod(size(V))))
end

fid     = fopen('gom_H_199601.bin','rb','ieee-be');
[SSHics,count]  = fread(fid,prod(size(SSH)),'float32');
fclose(fid);
if (count ~= prod(size(SSH)))
   error(sprintf('SSH record length wrong %d %d',count,prod(size(SSH))))
end

datvec = [Sics; Tics; Uics; Vics; SSHics];
disp(sprintf('total model size is %d',length(datvec)*8))

fid     = fopen('perfect_ics','wb','ieee-be');
fwrite(fid,  trec1,'int32');
fwrite(fid,seconds,'int32');
fwrite(fid,   days,'int32');
fwrite(fid,  trecN,'int32');

fwrite(fid,   rec1,'int32');
fwrite(fid, datvec,'float64');
fwrite(fid,   recN,'int32');
fclose(fid);

fid     = fopen('perfect_ics.txt','wt');
fprintf(fid,'%d %d\n',seconds,days);
fprintf(fid,'%.15e\n',datvec);
fclose(fid);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
