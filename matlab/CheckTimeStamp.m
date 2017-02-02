  function CheckTimeStamp(directory, varargin)
% Prints the timestamps of BINARY DART-format files
% and has option to print the advance-to-time (user must know if it exists)
%
% Inputs:
%   directory (string, required)  - directory containing .ics files
%   endian    (string, optional)  - 'ieee-be' ('big') or 'ieee-le' ('little'),
%                                   default is the native endian your system uses
%   fbase     (string, optional)  - base of ics files, default is 'filter_ics'
%   twoTimes  (string, optional)  - print advance-to time ('true'/'false'), default is 'false'
%
% If the dates are nonsensical, you probably have the wrong 'endian' specified.
%
% Outputs:
%   prints to screen filename and timestamp(s)
%
% Example 1: check all instances of filter_ics.nnnn using native endian-ness
%
%   CheckTimeStamp('/glade/scratch/syha/work_hires')
%
% Example 2: little-endian binary file 
%
%   directory = '/glade/scratch/syha/work_hires';
%   endian    = 'little';
%   CheckTimeStamp(directory,'endian', endian)
%
% Example 3: check all instances of assim_model_state_ic.nnnn for
%            both time stamps using native endianness
%
%   CheckTimeStamp('/glade/scratch/syha/work_hires', 'endian', 'native', ...
%                  'fbase', 'assim_model_state_ic', 'twoTimes', 'true')

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

timebase = datenum(1601,1,1); % start of Gregorian Calendar

p = inputParser; % create new instance of parser class
p.FunctionName = 'input parser :: requires input directory (string); valid optional inputs are, endian (string), fbase (string)';

% set defaults for optional parameters
defaultEndian    = 'native';
defaultFbase     = 'filter_ics';
defaulttwoTimes  = 'false';

addRequired(  p, 'directory', @ischar); % require directory, check input is a character array
addParamValue(p,    'endian', defaultEndian,   @ischar);
addParamValue(p,     'fbase', defaultFbase,    @ischar);
addParamValue(p,  'twoTimes', defaulttwoTimes, @ischar);

p.parse(directory, varargin{:}) % parse inputs

%% collect the results of parsing (makes code easier to read)

endianIn = p.Results.endian;
fbase    = p.Results.fbase;
twoTimes = p.Results.twoTimes;

% check inputs

assert(exist(directory, 'dir')==7, 'directory %s does not exist', directory)

tempstr  = strcat(fbase, '*');
ens_size = length(dir(fullfile(directory, tempstr)));

if (ens_size > 0)
   fprintf('The ensemble size is %d\n',ens_size)
else
   error('no %s files exist in %s',fbase,directory);
end

switch lower(endianIn)
   case {'big','ieee-be'}
      endian = 'ieee-be';
   case {'little','ieee-le'}
      endian = 'ieee-le';
   otherwise
      endian = 'native';
end

%% loop over each file - read the first 8 integers
%    If the file has only one timestamp record, the last
%    4 integers are garbage. rec1(1), rec(4), rec(5), and
%    rec(8) are record information only. Not useful.
%
%    Fortran record 1 ... [reclen][seconds][days][reclen]   (these are all int32)
%    Fortran record 2 ... [reclen][seconds][days][reclen]   (these are all int32)
%    Fortran record 3 ... [reclen][ model_state ][reclen]   (int32, N*real*8, int32)
%
%    Some files only have two records, in which case the situation is:
%
%    Fortran record 1 ... [reclen][seconds][days][reclen]   (these are all int32)
%    Fortran record 2 ... [reclen][ model_state ][reclen]   (int32, N*real*8, int32)
%
%    The time record closest to the model state is ALWAYS the valid_time of the model.

for i = 1:ens_size

   fname = sprintf('%s/%s.%04d',directory,fbase,i);

   if exist(fname, 'file')

      fid  = fopen(fname,'rb',endian);
      rec1 = fread(fid,8,'int32');
      fclose(fid);

      days    = rec1(3);
      seconds = rec1(2);
      fdays   = seconds/86400;
      advanceToTime = datestr(timebase + days + fdays);

      if strcmpi(twoTimes,'true')

         days2     = rec1(7);
         seconds2  = rec1(6);
         fdays2    = seconds2/86400;
         modelTime = datestr(timebase + days2 + fdays2);

         fprintf('%s.%04d advance-to-time %d %d (%s) & model time of %d %d (%s)\n'...
             ,fbase,i,days,seconds,advanceToTime, days2, seconds2, modelTime)

         % catch some, but not all errors
         if seconds2 < 0 || seconds2 > 86400
            error('%s might not have two timestamps', fname)
         end

      else

         fprintf('%s.%04d has timestamp of %d %d (%s)\n',fbase,i,days,seconds,advanceToTime)

      end

      if days < 0 || seconds < 0
         error('%s is incorrect endian, gives negative timestamp', endian)
      end

   else
      fprintf('WARNING : %s.%04d does not exist \n',fbase,i)
   end

end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
