function a = ReadBinaryObsSeq(fname,machineformat)
%% ReadBinaryObsSeq       reads the diagnostic output observation sequence file.
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

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin <= 1 )
   fname = 'obs_seq.final';
   machineformat = 'native';
elseif (nargin < 2 )
   fname = 'obs_seq.final';
end

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

fid = fopen(fname,'rb',machineformat);

% Read and parse the information that will help us read the rest of the file.

% Declare some integer types that will be used to match observations to
% types that have additional metadata and require special processing.

DOPPLER_RADIAL_VELOCITY = 12;        % default value for both formats
RAW_STATE_1D_INTEGRAL   = -99; oned_integral_already_read = 0;
GPSRO_REFRACTIVITY      = -99;


bogus1 = fread(fid,1,'int32');       % Read the first record header
values = char(fread(fid,[1 bogus1],'char'));  % Read the first record
bogus1 = fread(fid,1,'int32');       % Finish reading the first record


switch lower(values)
   case 'obs_sequence'      % pre-I or newer format

      bogus1 = fread(fid,1,'int32');       % Read the first record of write_obs_kind
      string = deblank(char(fread(fid,[1 bogus1],'char')));   % Read the first record
      bogus1 = fread(fid,1,'int32');       % Finish reading the first record

      bogus1  = fread(fid,1,'int32');
      numdefs = fread(fid,1,'int32')
      bogus1  = fread(fid,1,'int32');

      obskindnumber = zeros(1,numdefs);
      obskindstring = cell(1,numdefs);

      for idef =1:numdefs,

         bogus1  = fread(fid,1,'int32');   % record length
         tim     = fread(fid,1,'int32');

         string  = fread(fid,bogus1-2-2,'char');   % Read the character name
         tom     = deblank(char(string'));
         bogus1  = fread(fid,1,'int32');

         obskindnumber(idef) = tim;
         obskindstring(idef) = {tom};

         switch obskindstring{idef}
            case 'DOPPLER_RADIAL_VELOCITY'
               DOPPLER_RADIAL_VELOCITY = obskindnumber(idef);
            case 'RAW_STATE_1D_INTEGRAL'
               RAW_STATE_1D_INTEGRAL = obskindnumber(idef);
            case 'GPSRO_REFRACTIVITY'
               GPSRO_REFRACTIVITY = obskindnumber(idef);
         end
      end

      a = struct('filename',fname,'obskindnumber',obskindnumber);
      a.obskindstring = obskindstring;

otherwise    % before pre-I

   frewind(fid);
   a = struct('filename',fname);

end

%write(file_id) seq%num_copies, seq%num_qc, seq%num_obs, seq%max_num_obs

bogus1       = fread(fid,1,'int32');
num_copies   = fread(fid,1,'int32');
num_qc       = fread(fid,1,'int32');
num_obs      = fread(fid,1,'int32');
max_num_obs  = fread(fid,1,'int32');
bogus2       = fread(fid,1,'int32');

fprintf('num_copies  is %d\n',num_copies)
fprintf('num_qc      is %d\n',num_qc)
fprintf('num_obs     is %d\n',num_obs)
fprintf('max_num_obs is %d\n',max_num_obs)

% Read the metadata
for i = 1:num_copies,
   reclen = fread(fid,1,'int32');
   metadata{i} = deblank(char(fread(fid,[1 reclen],'char')));
   recend = fread(fid,1,'int32');
   if ( recend ~= reclen )
      str = sprintf('copy metadata read for copy %d failed %d %d', ...
                    i,reclen,recend);
      error(str)
   else
      fprintf('%s%s\n',metadata{i},'[END]')
   end
end

% Read the QC info
for i = 1:num_qc,

   reclen = fread(fid,1,'int32');
   qcdata{i} = deblank(char(fread(fid,[1 reclen],'char')));
   recend = fread(fid,1,'int32');

   if ( recend ~= reclen )
      str = sprintf('qc metadata read for copy %d failed %d %d', ...
                    i,reclen,recend);
      error(str)
   else
      fprintf('%s%s\n',qcdata{i},'[END]')
   end
end

% Read first/last time
bogus1 = fread(fid,1,'int32');
first_time = fread(fid,1,'int32');    % the indices of the adjacent obs
 last_time = fread(fid,1,'int32');    % the indices of the adjacent obs
bogus2 = fread(fid,1,'int32');

fprintf('first_time = %d  last_time = %d\n',first_time, last_time)

a.num_copies  = num_copies;
a.num_qc      = num_qc;
a.num_obs     = num_obs;
a.max_num_obs = max_num_obs;
a.first_time  = first_time;
a.last_time   = last_time;

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
radloc     = zeros(num_obs,3);
radwhich_vert = zeros(num_obs,1);
dir3d      = zeros(num_obs,3);

days(:)       = NaN;
secs(:)       = NaN;
evar(:)       = NaN;
kind(:)       = NaN;
prev_time(:)  = NaN;
next_time(:)  = NaN;
cov_group(:)  = NaN;
obs(:)        = NaN;
qc(:)         = NaN;
loc(:)        = NaN;
which_vert(:) = NaN;
radloc(:)     = NaN;
radwhich_vert(:) = NaN;
dir3d(:)      = NaN;

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
      fprintf('obs %d prev_time = %d next_time = %d cov_group = %d\n', ...
                    i,     prev_time(i),  next_time(i),  cov_group(i))
   end

   % Read the observation definition obs_def_mod:write_obs_def

   [loc(i,:) which_vert(i)] = read_location(fid,i);

   % Read the observation kind obs_def_mod:write_obs_def:write_kind

   bogus1  = fread(fid,1,'int32');
   kind(i) = fread(fid,1,'int32');   % obs_kind_type%index [integer]
   bogus2  = fread(fid,1,'int32');

   %------------------------------------------------------------
   % Specific kinds of observations have additional metadata.
   % Basically, we have to troll through the obs_def* files to
   % see what types need special attention.
   %------------------------------------------------------------

   if( kind(i) == DOPPLER_RADIAL_VELOCITY )

      [radloc(i,:) radwhich_vert(i)] = read_location(fid,i);

      dir3d(i,:) = read_orientation(fid,i);

      bogus1  = fread(fid,1,'int32');
      key     = fread(fid,1,'int32');
      bogus2  = fread(fid,1,'int32');

   elseif( kind(i) == RAW_STATE_1D_INTEGRAL )

      % Right now, I am throwing away the results.
      if ( ~ oned_integral_already_read )

         % the number of 1d_integral obs descriptions
         bogus1  = fread(fid,1,'int32');
         num_1d_integral_obs = fread(fid,1,'int32');
         bogus2  = fread(fid,1,'int32');

         % half_width, num_points, and localization_type for each
         for i1d=1:num_1d_integral_obs,
            bogus1  = fread(fid,1,'int32');
            half_width         = fread(fid,1,'float64'); % real(r8)
            num_points         = fread(fid,1,'int32');   % integer
            localization_type  = fread(fid,1,'int32');   % integer
            bogus2  = fread(fid,1,'int32');
            if (bogus1 ~= bogus2)
               error('RAW_STATE_1D_INTEGRAL: off track at %d %d',i,i1d)
            end
         end

         oned_integral_already_read = 1;  % true

      end

      bogus1  = fread(fid,1,'int32');
      key     = fread(fid,1,'int32');
      bogus2  = fread(fid,1,'int32');

   elseif( kind(i) == GPSRO_REFRACTIVITY )

      nada = read_gpsro_ref(fid,i);

   end

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

a.days       = days;
a.secs       = secs;
a.evar       = evar;
a.obs        = obs;
a.qc         = qc;
a.prev_time  = prev_time;
a.next_time  = next_time;
a.cov_group  = cov_group;
a.loc        = loc;
a.which_vert = which_vert;
a.kind       = kind;
a.metadata   = metadata;
a.qcdata     = qcdata;

fclose(fid);

%-------------------------------------------------------------------------------
% The companion helper functions ... trying to mimic the Fortran counterparts.
%-------------------------------------------------------------------------------

function [locs, which_vert] = read_location(fid,i)
% Read what was written by ANY location_mod:write_location
% Different location mods have different numbers of components to read.
% Must query the Fortran record information ...
%
% oned                 write(ifile)loc%x
% twod_sphere          write(ifile)loc%lon, loc%lat
% simple_threed_sphere write(ifile)loc%lon, loc%lat, loc%lev
% threed_sphere        write(ifile)loc%lon, loc%lat, loc%vloc, loc%which_vert
%
% Sort of cheating by returning the vectorized counterpart.
% Should do something with the determined dimension.

which_vert = 0;           % unless found to be otherwise
locs       = zeros(1,3);

reclen   = fread(fid,1,'int32');
if ( reclen == 8 )                          % oned
      locs(1) = fread(fid,1,'float64');
elseif ( reclen == 16 )                     % twod
      locs(1:2) = fread(fid,[1 2],'float64');
elseif ( reclen == 24 )                     % simple_threed_sphere
      locs(1:3) = fread(fid,[1 3],'float64');
elseif ( reclen == 28 )                     % threed_sphere
      locs(1:3) = fread(fid,[1 3],'float64');
      which_vert = fread(fid,1,'int32');
else
      error('Unknown record size for read_loc - got %d',reclen)
end
recend  = fread(fid,1,'int32');

if ( reclen ~= recend )
      error('off track at %d read_location %d ~= %d',i,reclen,recend)
end




function locs = read_orientation(fid,i)
% Read radar orientation

reclen  = fread(fid,1,'int32');
locs    = fread(fid,3,'float64');
recend  = fread(fid,1,'int32');

if ( reclen ~= recend )
      error('off track at %d read_orientation %d ~= %d',i,reclen,recend)
end




function key = read_gpsro_ref(fid,i)
% obs_def_gps_mod:read_gpsro_ref    equivalent.
%
% for now, no return values.
%
% if it dies in here, it dies everywhere, so we don't need a return code.
%
%   read(ifile) key
%   read(ifile) rfict(key), step_size(key), ray_top(key), &
%              (ray_direction(ii, key), ii=1, 3), gpsro_ref_form(key)
%
% real(r8) :: rfict, step_size, ray_top, ray_direction
% character(len=6) :: gpsro_ref_form

reclen  = fread(fid,1,'int32');
key     = fread(fid,1,'int32');
recend  = fread(fid,1,'int32');

reclen  = fread(fid,1,'int32');
rfict          = fread(fid,1,'float64');
step_size      = fread(fid,1,'float64');
ray_top        = fread(fid,1,'float64');
ray_direction  = fread(fid,3,'float64');
string         = fread(fid,6,'char');
recend  = fread(fid,1,'int32');

gpsro_ref_form = char(string');  % convert binary to ascii string


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
