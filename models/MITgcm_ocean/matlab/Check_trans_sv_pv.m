function Check_trans_sv_pv
%% Check_trans_sv_pv()  verify the packing/unpacking of the DART state vector
% there are a lot of hardwired values ...

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%time for timestep 0 day=128565, sec=21600
%time for timestep 0 1953 Jan 01 06:00:00
%time for 0000001344 day=128992, sec=21600
%time for 0000001344 1954 Mar 04 06:00:00

dt = 900;
startDate_1=19960101;
startDate_2=60000;
tbase = datenum(1996,1,1,6,0,0);

S   = rdmds(  '~/MITgcm/data40levels/S.0000001344');
T   = rdmds(  '~/MITgcm/data40levels/T.0000001344');
U   = rdmds(  '~/MITgcm/data40levels/U.0000001344');
V   = rdmds(  '~/MITgcm/data40levels/V.0000001344');
SSH = rdmds('~/MITgcm/data40levels/Eta.0000001344');

modelsize = prod(size(S)) + prod(size(T)) + ...
            prod(size(U)) + prod(size(V)) + prod(size(SSH));

toffset = (01344*dt)/86400

[nx ny nz] = size(S);

fid = fopen('assim_model_state_ud','rb','ieee-be');

rec1    = fread(fid,1,'int32');
seconds = fread(fid,1,'int32');
days    = fread(fid,1,'int32');
recN    = fread(fid,1,'int32');

if (rec1 ~= recN) 
   error('first record mismatch')
end
disp(sprintf('time tag is days/seconds : %d %d',days,seconds))

recN    = fread(fid,        1,'int32');
datmat  = fread(fid,modelsize,'float64');
recN    = fread(fid,        1,'int32');
fclose(fid);

fid = fopen('assim_model_state_ic','wb','ieee-be');
fwrite(fid,rec1,'int32');
fwrite(fid,seconds+900,'int32');
fwrite(fid,days+2,'int32');
fwrite(fid,rec1,'int32');

fwrite(fid,rec1,'int32');
fwrite(fid,seconds,'int32');
fwrite(fid,days,'int32');
fwrite(fid,rec1,'int32');

fwrite(fid,  recN,'int32');
fwrite(fid,datmat,'float64');
fwrite(fid,  recN,'int32');
fclose(fid);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
