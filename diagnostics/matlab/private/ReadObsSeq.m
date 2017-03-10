function a = ReadObsSeq(fname)
%% ReadObsSeq       reads the diagnostic output observation sequence file.
%
% The observation sequence file can be ascii or binary -- and either-endian
% if they are binary. A couple quick checks figure out the file format.
%
% a = ReadObsSeq('obs_seq.final');
%
% there are many returned components of 'a' ...
% >> fieldnames(a)
%ans =
%   'filename'
%   'num_copies'
%   'num_qc'
%   'num_obs'
%   'max_num_obs'
%   'first_time'
%   'last_time'
%   'days'
%   'secs'
%   'evar'
%   'obs'
%   'qc'
%   'prev_time'
%   'next_time'
%   'cov_group'
%   'loc'
%   'which_vert'
%   'kind'
%
% Uses ReadASCIIObsSeq.m  and ReadBINARYObsSeq.m

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin < 1 )
   fname = 'obs_seq.final';
end

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

% Determine if the file is an ascii file:

fid    = fopen(fname,'rt');
aline  = fgetl(fid);
values = sscanf(aline,'%s');
fclose(fid);

if (strcmp(values,'obs_sequence'))
   flavor = 'ASCII';
end

% Determine if the file is a little-endian binary file:

fid    = fopen(fname,'rb','ieee-le');
bogus1 = fread(fid,1,'int32');
fclose(fid);

if ( bogus1 == 4*4 )
   flavor = 'ieee-le';
end

% Determine if the file is a big-endian binary file:

fid    = fopen(fname,'rb','ieee-be');
bogus1 = fread(fid,1,'int32');
fclose(fid);

if ( bogus1 == 4*4 )
   flavor = 'ieee-be';
end

% presumably, now we know ...

switch  lower(flavor)
   case 'ascii'
	   a = ReadASCIIObsSeq(fname);
   case 'ieee-le'
	   a = ReadBinaryObsSeq(fname,'ieee-le');
   case 'ieee-be'
	   a = ReadBinaryObsSeq(fname,'ieee-be');
   otherwise
      error('Unable to determine format of %s',fname)
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
