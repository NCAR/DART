function a = ReadASCIIObsSeq(fname)
% ReadASCIIObsSeq       reads the diagnostic output observation sequence file.
%
% 
% a = ReadASCIIObsSeq('obs_seq.final');
%
% This is pretty slow -- lots of logic and nested loops -- hard to vectorize. 

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

%
% Intel F90 has a nasty habit of writing numbers like: 5.44307648074716D-003
% Matlab does not interpret the D correctly, it likes 'E' ... so all
% the scans for reals have to convert any possible 'D's to 'E's ... slowness.
%

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

% Read the metadata -- just throwing it away for now
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

   aline      = fgetl(fid);
   values     = sscanf(aline,'%*s %d');
   key        = values(1);

   if ( key    ~= i ) 
      error(sprintf('key %d is not %d -- and it should be ...',i,key))
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
   obdefstr = sscanf(aline,'%s');

   if ( ~ strcmp(obdefstr,'obdef') )  
      error(sprintf(' read %s instead of ''obdef''',obdefstr))
   end

   % Read what was written by location_mod:write_location

   [loc(i,:) which_vert(i)] = read_location(fid,i); 
   kind(i)            = read_kind(fid,i);  % companion to obs_kind_mod:write_kind
   [days(i), secs(i)] = read_time(fid);    % companion to time_manager_mod:write_time

   % And finally: the error variance

   aline    = cvrtline(fgetl(fid)); 
   evar(i)  = sscanf(aline,'%e');

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
   kindstr = sscanf(aline,'%s');
   if ( ~ strcmp(kindstr,'kind') )  
      error(sprintf(' read %s instead of ''kind'' at obs %d',obdefstr,i))
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
         aline    = fgetl(fid);
         values   = sscanf(aline,'%e');
         locs(1)  = values(1);
      case 'loc2s'       % location/twod_sphere/location_mod.f90:write_location
         aline    = fgetl(fid);
         values   = sscanf(aline,'%e');
	 %lon      = values(1);
	 %lat      = values(2);
	 %lev      = values(3);
         locs(1:2) = values(1:2);
      case 'loc3s'       % location/simple_threed_sphere/location_mod.f90:write_location
         aline    = fgetl(fid);
         values   = sscanf(aline,'%e');
	 %lon      = values(1);
	 %lat      = values(2);
	 %lev      = values(3);
         locs(1:3) = values(1:3);
      case 'loc3d'       % location/threed_sphere/location_mod.f90:write_location
         aline    = fgetl(fid);
         values   = sscanf(aline,'%e');
         locs(1:3) = values(1:3);
	 %lon        = values(1);
	 %lat        = values(2);
	 %lev        = values(3);
	 which_vert = values(4);  % which kind of vertical coord system being used.
      otherwise
         error(sprintf('unrecognized location type ... read %s',locstr))
end
