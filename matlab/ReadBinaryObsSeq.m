function a = ReadBinaryObsSeq(fname,machineformat)
% ReadBinaryObsSeq       reads the diagnostic output observation sequence file.
% 
% machineformat = 'n';
% a = ReadBinaryObsSeq('obs_seq.final',machineformat)
%
%    MACHINEFORMAT is one of the following strings:
% 
%    'native'      or 'n' - local machine format - the default
%    'ieee-le'     or 'l' - IEEE floating point with little-endian
%                           byte ordering
%    'ieee-be'     or 'b' - IEEE floating point with big-endian
%                           byte ordering
%    'vaxd'        or 'd' - VAX D floating point and VAX ordering
%    'vaxg'        or 'g' - VAX G floating point and VAX ordering
%    'cray'        or 'c' - Cray floating point with big-endian
%                           byte ordering
%    'ieee-le.l64' or 'a' - IEEE floating point with little-endian
%                           byte ordering and 64 bit long data type
%    'ieee-be.l64' or 's' - IEEE floating point with big-endian byte
%                           ordering and 64 bit long data type.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

if (nargin <= 1 )
   fname = 'obs_seq.final';
   machineformat = 'native';
elseif (nargin < 2 )
   fname = 'obs_seq.final';
end

fid = fopen(fname,'rb',machineformat);

% Read and parse the information that will help us read the rest of the file.

%write(file_id) seq%num_copies, seq%num_qc, seq%num_obs, seq%max_num_obs

bogus1       = fread(fid,1,'int32');
num_copies   = fread(fid,1,'int32');
num_qc       = fread(fid,1,'int32');
num_obs      = fread(fid,1,'int32');
max_num_obs  = fread(fid,1,'int32');
bogus2       = fread(fid,1,'int32');

disp(sprintf('num_copies  is %d',num_copies))
disp(sprintf('num_qc      is %d',num_qc))
disp(sprintf('num_obs     is %d',num_obs))
disp(sprintf('max_num_obs is %d',max_num_obs))

% Read the metadata -- just throwing it away for now.
for i = 1:num_copies,

   reclen = fread(fid,1,'int32');
   metadatastr = char(fread(fid,[1 reclen],'char'));
   recend = fread(fid,1,'int32');

   if ( recend ~= reclen ) 
      str = sprintf('copy metadata read for copy %d failed %d %d', ...
                    i,reclen,recend);
      error(str)
   else
      disp(sprintf('%s%s',metadatastr,'[END]'))
   end
end

% Read the QC info
for i = 1:num_qc,

   reclen = fread(fid,1,'int32');
   qcdata = char(fread(fid,[1 reclen],'char'));
   recend = fread(fid,1,'int32');

   if ( recend ~= reclen ) 
      str = sprintf('qc metadata read for copy %d failed %d %d', ...
                    i,reclen,recend);
      error(str)
   else
      disp(sprintf('%s%s',qcdata,'[END]'))
   end
end

% Read first/last time 
bogus1 = fread(fid,1,'int32');
first_time = fread(fid,1,'int32');    % the indices of the adjacent obs
 last_time = fread(fid,1,'int32');    % the indices of the adjacent obs
bogus2 = fread(fid,1,'int32');

disp(sprintf('first_time = %d  last_time = %d',first_time, last_time))

% Read the observations
days       = zeros(num_obs,1);
secs       = zeros(num_obs,1);
evar       = zeros(num_obs,1);
kind       = zeros(num_obs,1);
prev_time  = zeros(num_obs,1);
next_time  = zeros(num_obs,1);
cov_group  = zeros(num_obs,1);
obs        = zeros(num_copies,num_obs);
qc         = zeros(num_qc,num_obs);
loc        = zeros(num_obs,3);  % declare for worst case, pare down later.
which_vert = zeros(num_obs,1);

for i = 1:num_obs, % Read what was written by obs_sequence_mod:write_obs

   for j = 1:num_copies,
      bogus1   = fread(fid,1,'int32');
      obs(j,i) = fread(fid,1,'float64');
      bogus2   = fread(fid,1,'int32');
   end

   for j = 1:num_qc,
      bogus1   = fread(fid,1,'int32');
      qc(j,i)  = fread(fid,1,'float64');
      bogus2   = fread(fid,1,'int32');
   end
    
   bogus1   = fread(fid,1,'int32');
   prev_time(i) = fread(fid,1,'int32');   % obs_type%prev_time
   next_time(i) = fread(fid,1,'int32');   % obs_type%next_time
   cov_group(i) = fread(fid,1,'int32');   % obs_type%cov_group
   bogus2   = fread(fid,1,'int32');

   if ( mod(i,1000) == 0 ) 
      disp(sprintf('obs %d prev_time = %d next_time = %d cov_group = %d', ...
                    i,     prev_time(i),  next_time(i),  cov_group(i)))
   end

   % Read the observation definition obs_def_mod:write_obs_def

   % Read what was written by location_mod:write_location
   % Here is the tricky part. Different location mods have 
   % different numbers of components to read here. Must query
   % the Fortran record information ...
   % oned                 write(ifile)loc%x
   % twod_sphere          write(ifile)loc%lon, loc%lat
   % simple_threed_sphere write(ifile)loc%lon, loc%lat, loc%lev
   % threed_sphere        write(ifile)loc%lon, loc%lat, loc%vloc, loc%which_vert

   reclen   = fread(fid,1,'int32');
   if ( reclen == 8 ) 
      loc(i,  1) = fread(fid,1,'float64');
   elseif ( reclen == 16 ) 
      loc(i,1:2) = fread(fid,[1 2],'float64');
   elseif ( reclen == 24 ) 
      loc(i,1:3) = fread(fid,[1 3],'float64');
   elseif ( reclen == 28 ) 
      loc(i,1:3) = fread(fid,[1 3],'float64');
      which_vert(i) = fread(fid,1,'int32');
   else
      error(sprintf('Unknown record size for read_loc - got %d',reclen))
   end
   recend  = fread(fid,1,'int32');

   if ( reclen ~= recend ) 
      error(sprintf('off track after read_location %d ~= %d',reclen,recend))
   end

   % Read the observation kind obs_def_mod:write_obs_def:write_kind

   bogus1  = fread(fid,1,'int32');
   kind(i) = fread(fid,1,'int32');   % obs_kind_type%index [integer]
   bogus2  = fread(fid,1,'int32');

   % Read the time                  time_manager:write_time

   bogus1  = fread(fid,1,'int32');
   secs(i) = fread(fid,1,'int32');
   days(i) = fread(fid,1,'int32');
   bogus2  = fread(fid,1,'int32');

   % And finally: the error variance

   bogus1  = fread(fid,1,'int32');
   evar(i) = fread(fid,1,'float64');   % obs_kind_type%index [integer]
   bogus2  = fread(fid,1,'int32');

end

a = struct( 'filename'   , fname       , ...
	    'num_copies' , num_copies  , ...
            'num_qc'     , num_qc      , ...
	    'num_obs'    , num_obs     , ...
	    'max_num_obs', max_num_obs , ...
	    'first_time' , first_time  , ...
	    'last_time'  , last_time   , ...
	    'days'       , days        , ...
	    'secs'       , secs        , ...
	    'evar'       , evar        , ...
	    'obs'        , obs         , ...
	    'qc'         , qc          , ...
	    'prev_time'  , prev_time   , ...
	    'next_time'  , next_time   , ...
	    'cov_group'  , cov_group   , ...
	    'loc'        , loc         , ...
	    'which_vert' , which_vert  , ...
	    'kind'       , kind        );
%    'metadata', metadata, ...
%    'qcdata', qcdata, ...

fclose(fid);

%-------------------------------------------------------------------------------
% The companion helper functions ... trying to mimic the Fortran counterparts.
%-------------------------------------------------------------------------------
