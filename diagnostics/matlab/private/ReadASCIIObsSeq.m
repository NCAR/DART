function a = ReadASCIIObsSeq(fname)
%% ReadASCIIObsSeq       reads the diagnostic output observation sequence file.
%
%
% a = ReadASCIIObsSeq('obs_seq.final');
%
% This is pretty slow -- lots of logic and nested loops -- hard to vectorize.
% The best thing is probably to turn it into a '.mex' file.
%
% The file contains a linked list which we are reading sequentially.
% The resulting sequence "a.obs" is NOT guaranteed to be in a temporally ascending order.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin < 1 )
   fname = 'obs_seq.final';
end

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

fid = fopen(fname,'rt');

%
% Intel F90 has a nasty habit of writing numbers like: 5.44307648074716D-003
% Matlab does not interpret the D correctly, it likes 'E' ... so all
% the scans for reals have to convert any possible 'D's to 'E's ... slowness.
%

% Read and parse the information that will help us read the rest of the file.

% First, determine if it is 'old school' (i.e. pre-I)
% If the first line contains the word 'obs_sequence', it is pre-I or better.
DOPPLER_RADIAL_VELOCITY = 12;                      % default value for both formats
RAW_STATE_1D_INTEGRAL   = -99; oned_integral_already_read = 0;
GPSRO_REFRACTIVITY      = -99;

aline      = fgetl(fid);
values     = sscanf(aline,'%s %*');   % Just concerned with the first word.

switch lower(values)
   case 'obs_sequence'      % pre-I or newer format

      aline = fgetl(fid);   % skip line containing the word 'obs_kind_definitions'
      aline = fgetl(fid);
      numdefs = sscanf(aline,'%d');

      obskindnumber = zeros(1,numdefs);
      obskindstring = cell(1,numdefs);

      for idef =1:numdefs,
         aline = fgetl(fid);

         tim = sscanf(aline,'%d %*s');
         tom = deblank(sscanf(aline,'%*s %s'));

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

      aline      = fgetl(fid);
      values     = sscanf(aline,'%*s %d %*s %d');
      num_copies = values(1);
      num_qc     = values(2);

      a = struct('filename',fname,'obskindnumber',obskindnumber);
      a.obskindstring = obskindstring;

   otherwise                % before pre-I

      values     = sscanf(aline,'%*s %d %*s %d');
      num_copies = values(1);
      num_qc     = values(2);

      a = struct('filename',fname);

end

aline        = fgetl(fid);
values       = sscanf(aline,'%*s %d %*s %d');
num_obs      = values(1);
max_num_obs  = values(2);

disp(sprintf('num_copies  is %d',num_copies))
disp(sprintf('num_qc      is %d',num_qc))
disp(sprintf('num_obs     is %d',num_obs))
disp(sprintf('max_num_obs is %d',max_num_obs))

% Read the metadata
for i = 1:num_copies,
   metadata{i} = deblank(fgetl(fid));
end

% Read the QC info
for i = 1:num_qc,
   qcdata{i} = deblank(fgetl(fid));
end

% Read and parse the information that will help us read the rest of the file.
aline      = fgetl(fid);
values     = sscanf(aline,'%*s %d %*s %d');
first_time = values(1);
last_time  = values(2);

disp(sprintf('first_time = %d  last_time = %d',first_time, last_time))

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

days(:)          = NaN;
secs(:)          = NaN;
evar(:)          = NaN;
kind(:)          = NaN;
prev_time(:)     = NaN;
next_time(:)     = NaN;
cov_group(:)     = NaN;
obs(:)           = NaN;
qc(:)            = NaN;
loc(:)           = NaN;
which_vert(:)    = NaN;
radloc(:)        = NaN;
radwhich_vert(:) = NaN;
dir3d(:)         = NaN;

for i = 1:num_obs, % Read what was written by obs_sequence_mod:write_obs

   aline      = fgetl(fid);
   values     = sscanf(aline,'%*s %d');
   key        = values(1);

   if ( key    ~= i )
      error(sprintf('key %d is not %d -- and it should be ...',i,key))
   % else
   %    disp(sprintf('reading observation %d',i))
   end

   for j = 1:num_copies,
      aline    = cvrtline(fgetl(fid));
      obs(j,i) = sscanf(aline,'%e');
   end

   for j = 1:num_qc,
      aline    = cvrtline(fgetl(fid));
      qc(j,i)  = sscanf(aline,'%e');
   end

   aline     = fgetl(fid);
   values    = sscanf(aline,'%d');
   prev_time(i) = values(1);
   next_time(i) = values(2);
   cov_group(i) = values(3);

   if ( mod(i,1000) == 0 )
      disp(sprintf('obs %d prev_time = %d next_time = %d cov_group = %d', ...
                    i,     prev_time(i),  next_time(i),  cov_group(i)))
   end

   % Read the observation definition obs_def_mod:write_obs_def

   aline    = fgetl(fid);
   obdefstr = deblank(sscanf(aline,'%s'));

   if ( ~ strcmp(obdefstr,'obdef') )
      error(sprintf(' read %s instead of ''obdef''',obdefstr))
   end

   % Read what was written by location_mod:write_location

   [loc(i,:) which_vert(i)] = read_location(fid,i);
   kind(i)                  = read_kind(fid,i);  % companion to obs_kind_mod:write_kind

   %------------------------------------------------------------
   % Specific kinds of observations have additional metadata.
   % Basically, we have to troll through the obs_def* files to
   % see what types need special attention.
   %------------------------------------------------------------

   if(kind(i) == DOPPLER_RADIAL_VELOCITY)

      read_platform(fid,i);
      [radloc(i,:) radwhich_vert(i)] = read_location(fid,i);
      dir3d                          = read_orientation(fid,i);
      aline      = fgetl(fid);     % read the line with the key

   elseif(kind(i) == RAW_STATE_1D_INTEGRAL )

      % Right now, I am throwing away the results.
      if ( ~ oned_integral_already_read )

         aline = fgetl(fid);  % the number of 1d_integral obs descriptions
         num_1d_integral_obs = sscanf(aline,'%d');

         for i1d=1:num_1d_integral_obs,
            aline = fgetl(fid);  % half_width, num_points, and localization_type for each
         end

         oned_integral_already_read = 1;  % true

      end

      aline = fgetl(fid);  % get obs_def key

   elseif( kind(i) == GPSRO_REFRACTIVITY )

      read_gpsro_ref(fid,i);

   end

   [days(i), secs(i)] = read_time(fid);    % companion to time_manager_mod:write_time

   % And finally: the error variance

   aline    = cvrtline(fgetl(fid));
   evar(i)  = sscanf(aline,'%e');

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
if ( num_copies > 0); a.metadata = metadata; end;
if ( num_qc     > 0); a.qcdata   = qcdata; end;

fclose(fid);

%-------------------------------------------------------------------------------
% The companion helper functions ... trying to mimic the Fortran counterparts.
%-------------------------------------------------------------------------------

function ostr = cvrtline(istring)
% Try to replace all the stupid 'D's with 'E's
ostr = istring;
inds = find(ostr == 'D');
if ~ isempty(inds)
   ostr(inds) = 'E';
end


function [days, secs] = read_time(fid)
   aline  = fgetl(fid);
   values = sscanf(aline,'%d');
   secs   = values(1);
   days   = values(2);



function kind = read_kind(fid,i)
   aline   = fgetl(fid);
   kindstr = deblank(sscanf(aline,'%s'));
   if ( ~ strcmp(kindstr,'kind') )
      error(sprintf(' read %s instead of ''kind'' at obs %d',kindstr,i))
   end
   aline   = fgetl(fid);
   kind    = sscanf(aline,'%d');



function [locs, which_vert] = read_location(fid,i)
% Read what was written by ANY location_mod:write_location
%
% Sort of cheating by returning the vectorized counterpart.
% Should do something with the determined dimension.

which_vert = 0;           % unless found to be otherwise
locs       = zeros(1,3);

aline    = fgetl(fid);
locstr   = sscanf(aline,'%s');

switch lower(locstr)
      case 'loc1d'       % location/oned/location_mod.f90:write_location
         aline    = cvrtline(fgetl(fid));
         values   = sscanf(aline,'%e');
         locs(1)  = values(1);
      case 'loc2s'       % location/twod_sphere/location_mod.f90:write_location
         aline    = cvrtline(fgetl(fid));
         values   = sscanf(aline,'%e');
	 %lon      = values(1);
	 %lat      = values(2);
         locs(1:2) = values(1:2);
      case 'loc3s'       % location/simple_threed_sphere/location_mod.f90:write_location
         aline    = cvrtline(fgetl(fid));
         values   = sscanf(aline,'%e');
	 %lon      = values(1);
	 %lat      = values(2);
	 %lev      = values(3);
         locs(1:3) = values(1:3);
      case 'loc3d'       % location/threed_sphere/location_mod.f90:write_location
         aline    = cvrtline(fgetl(fid));
         values   = sscanf(aline,'%e');
         locs(1:3) = values(1:3);
	 %lon        = values(1);
	 %lat        = values(2);
	 %lev        = values(3);
	 which_vert = values(4);  % which kind of vertical coord system being used.
         %aline      = cvrtline(fgetl(fid));
         %which_vert = sscanf(aline,'%d');

      otherwise
         error(sprintf('unrecognized location type ... read %s',locstr))
end


function read_platform(fid,i)
% obs_def_radar_mod:write_rad_vel    component.

aline       = fgetl(fid);
platformstr = deblank(sscanf(aline,'%s'));
if ( ~ strcmp(platformstr,'platform') )
   error(sprintf(' read %s instead of ''platform'' at obs %d',platformstr,i))
end



function dir3d = read_orientation(fid,i)
% obs_def_radar_mod:write_rad_vel    component.

aline    = fgetl(fid);
platformstr = deblank(sscanf(aline,'%s'));
if ( ~ strcmp(platformstr,'dir3d') )
    error(sprintf(' read %s instead of ''dir3d''',platformstr))
end
aline    = cvrtline(fgetl(fid));
values   = sscanf(aline,'%e');
dir3d    = values(1:3);



function read_gpsro_ref(fid,i)
% obs_def_gps_mod:read_gpsro_ref    equivalent.
%
% for now, no return values.
%
% if it dies in here, it dies everywhere, so we don't need a return code.
%
% read(ifile, FMT='(a8, i8)') header, key
% read(ifile, *) rfict(key), step_size(key), ray_top(key), &
%                  (ray_direction(ii, key), ii=1, 3), gpsro_ref_form(key)
%
% real(r8) :: rfict, step_size, ray_top, ray_direction
% character(len=6) :: gpsro_ref_form

aline  = fgetl(fid);
gpskey = deblank(sscanf(aline,'%s %*d'));
if ( ~ strcmp(gpskey,'gpsroref') )
    error(sprintf(' read %s instead of ''gpsroref''',platformstr))
end
aline  = cvrtline(fgetl(fid));
values = sscanf(aline,'%e %*s');
form   = sscanf(aline,'%*e %s');
if (length(values) ~= 6)
   error(sprintf('read %d items instead of 6 for obs# %d',length(values),i))
end
rfict         = values(1);
step_size     = values(2);
ray_top       = values(3);
ray_direction = values(4:6);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
