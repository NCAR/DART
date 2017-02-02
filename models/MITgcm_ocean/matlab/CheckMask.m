function b = CheckMask()
%% CheckMask 
% 
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%mitbase = '/fs/image/home/nancy/subversion/trunk/models/MITgcm_ocean/data2/';

mitbase = '/ptmp/thoar/MITgcm/exp2';
fname = sprintf('%s/advance_temp0/bathymetry.bin',mitbase);
b     = rdslab(fname);

%Sinds = find(S == 0.0); clear S;
%
%T     = rdmds(sprintf('%s/%s.0000040992',mitbase,'T'));
%Tinds = find(T == 0.0); clear T;
%
%U     = rdmds(sprintf('%s/%s.0000040992',mitbase,'U'));
%Uinds = find(U == 0.0); clear U;
%
%V     = rdmds(sprintf('%s/%s.0000040992',mitbase,'V'));
%Vinds = find(V == 0.0); clear V;
%
%SSH = rdmds(sprintf('%s/%s.0000040992',mitbase,'Eta'));

%disp(sprintf('S has %d zeros',length(Sinds)))
%disp(sprintf('T has %d zeros',length(Tinds)))
%disp(sprintf('U has %d zeros',length(Uinds)))
%disp(sprintf('V has %d zeros',length(Vinds)))

function x = rdslab(fname)

if (exist(fname,'file'))
   disp(sprintf('Opening %s',fname))
else
   error('Opening %s',fname)
end

fid       = fopen(fname,'rb','ieee-be');
[x count] = fread(fid,[256 225],'float32');

if (count ~= 256*225)
   error('only read %d of %d items from %s',count,256*225,fname)
end
fclose(fid)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
