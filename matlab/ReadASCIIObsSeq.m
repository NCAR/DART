function a = ReadASCIIObsSeq(fname)
% ReadASCIIObsSeq       reads the diagnostic output observation sequence file.
%
% 
% a = ReadASCIIObsSeq('obs_seq.final')
%

% binary ones might have endian issues.

% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% $Source$
% $Revision$
% $Date$

if (nargin < 1 )
   fname = 'obs_seq.final';
end

fid = fopen(fname,'rt');

% Read and parse the information that will help us read the rest of the file.
aline      = fgetl(fid);
values     = sscanf(aline,'%*s %d %*s %d');
num_copies   = values(1);
num_qc       = values(2);

aline      = fgetl(fid);
values     = sscanf(aline,'%*s %d %*s %d');
num_obs      = values(1);
max_num_obs  = values(2);

disp(sprintf('num_copies  is %d',num_copies))
disp(sprintf('num_qc      is %d',num_qc))
disp(sprintf('num_obs     is %d',num_obs))
disp(sprintf('max_num_obs is %d',max_num_obs))

% Read the metadata
for i = 1:num_copies,
   metadata{i} = fgetl(fid);
end

% Read the QC info
for i = 1:num_qc,
   qcdata{i} = fgetl(fid);
end

% Read and parse the information that will help us read the rest of the file.
aline      = fgetl(fid);
values     = sscanf(aline,'%*s %d %*s %d');
first_time = values(1);
last_time  = values(2);

% Read the observations
days = zeros(num_obs,1);
secs = zeros(num_obs,1);
evar = zeros(num_obs,1);
kind = zeros(num_obs,1);
key  = zeros(num_obs,1);
obs  = zeros(num_copies,num_obs);
qc   = zeros(num_qc,num_obs);

%for i = 1:num_obs, % Read what was written by obs_sequence_mod:write_obs
for i = 1:10, % Read what was written by obs_sequence_mod:write_obs

   aline      = fgetl(fid);
   values     = sscanf(aline,'%*s %d');
   key(i)     = values(1);

   if ( key(i) ~= i ) 
      error(sprintf('key %d is not %d -- and it should be ...',i,key(i)))
   end
 
   for j = 1:num_copies,
      aline    = fgetl(fid);
      obs(j,i) = sscanf(aline,'%f');
   end

   for j = 1:num_qc,
      aline    = fgetl(fid);
      qc(j,i)  = sscanf(aline,'%f');
   end
    
   aline     = fgetl(fid);
   values    = sscanf(aline,'%d');
   prev_time = values(1);
   next_time = values(2);
   cov_group = values(3);

   % Read the observation definition obs_def_mod:write_obs_def

   aline    = fgetl(fid);
   obdefstr = sscanf(aline,'%s');

   if ( ~ strcmp(obdefstr,'obdef') )  
      error(sprintf(' read %s instead of ''obdef''',obdefstr))
   end

   % Read what was written by location_mod:write_location

   location = read_location(fid,i);

%   aline    = fgetl(fid);
%   locstr   = sscanf(aline,'%s'); 
%
%   switch lower(locstr)
%      case 'loc1d'       % location/oned/location_mod.f90:write_location
%         aline    = fgetl(fid);
%         location = sscanf(aline,'%f');
%      case 'loc2s'       % location/twod_sphere/location_mod.f90:write_location
%                         % and unfortunately simple_threed_sphere ...
%         aline    = fgetl(fid);
%         values   = sscanf(aline,'%f');
%	 lon      = values(1);
%	 lat      = values(2);
%	 if ( length(values) > 2 )
%	    lev      = values(3);
%	 end 
%      case 'loc3d'       % location/threed_sphere/location_mod.f90:write_location
%         aline    = fgetl(fid);
%         values   = sscanf(aline,'%f');
%	 lon        = values(1);
%	 lat        = values(2);
%	 lev        = values(3);
%	 which_vert = values(4);  % which kind of vertical coord system being used.
%      otherwise
%         error(sprintf('unrecognized location type ... read %s',locstr))
%   end

   kind(i)            = read_kind(fid,i);  % companion to obs_kind_mod:write_kind 
   [days(i), secs(i)] = read_time(fid);    % companion to time_manager_mod:write_time

   % And finally: the error variance

   aline    = fgetl(fid);
   evar(i)  = sscanf(aline,'%f');

end

a = struct( 'filename',fname, ...
            'num_qc',num_qc, ...
	    'num_copies', num_copies, ...
	    'num_obs',num_obs, ...
	    'max_num_obs', max_num_obs, ...
	    'first_time', first_time, ...
	    'last_time', last_time, ...
	    'key' , key,  ...
	    'days', days, ...
	    'secs', secs, ...
	    'evar', evar, ...
	    'obs',  obs,  ...
	    'qc',   qc,   ...
	    'kind', kind);
%    'metadata', metadata, ...
%    'qcdata', qcdata, ...

fclose(fid);

%-------------------------------------------------------------------------------
% The companion helper functions ... trying to mimic the Fortran counterparts.
%-------------------------------------------------------------------------------

function [days, secs] = read_time(fid)
   aline  = fgetl(fid);
   values = sscanf(aline,'%d');
   secs   = values(1);
   days   = values(2);



function kind = read_kind(fid,i)
   aline   = fgetl(fid);
   kindstr = sscanf(aline,'%s');
   if ( ~ strcmp(kindstr,'kind') )  
      error(sprintf(' read %s instead of ''kind'' at obs %d',obdefstr,i))
   end
   aline   = fgetl(fid);
   kind    = sscanf(aline,'%d');



function values = read_location(fid,i)
% Read what was written by ANY location_mod:write_location
%
% Sort of cheating by returning the vectorized counterpart.
% Should do something with the determined dimension.
aline    = fgetl(fid);
locstr   = sscanf(aline,'%s'); 

switch lower(locstr)
      case 'loc1d'       % location/oned/location_mod.f90:write_location
         aline    = fgetl(fid);
         values   = sscanf(aline,'%f');
         loc      = values(1);
      case 'loc2s'       % location/twod_sphere/location_mod.f90:write_location
                         % and unfortunately simple_threed_sphere ...
         aline    = fgetl(fid);
         values   = sscanf(aline,'%f');
	 lon      = values(1);
	 lat      = values(2);
	 if ( length(values) > 2 )
	    lev      = values(3);
	 end 
      case 'loc3d'       % location/threed_sphere/location_mod.f90:write_location
         aline    = fgetl(fid);
         values   = sscanf(aline,'%f');
	 lon        = values(1);
	 lat        = values(2);
	 lev        = values(3);
	 which_vert = values(4);  % which kind of vertical coord system being used.
      otherwise
         error(sprintf('unrecognized location type ... read %s',locstr))
end
