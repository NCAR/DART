function pinfo_out = CheckModelCompatibility(arg1, arg2);
% CheckModelCompatibility tries to ensure that two netcdf files can be compared.
% There are 2 ways to call this:  with 2 filenames, or with an already existing
% pinfo struct (with 2 filenames and 2 2-vector arrays for start/stop times).
% This routine fills in the 2-vectors with the time overlap region in a
% pinfo struct.
% If the time indices are common between the 2 files it returns [1,nlength]
% for both; if no overlap [-1,-1] for both; otherwise the [start,end] indices 
% for each array.

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$


if (nargin == 1)      % better be a pinfo struct with at least these fields
  file1 = arg1.truth_file;  % string
  file2 = arg1.diagn_file;  % string
  time1 = arg1.truth_time;  % [a,b] array
  time2 = arg1.diagn_time;  % [a,b] array
  pinfo_out = arg1;
elseif (nargin == 2)  % pair of filenames
  file1 = arg1;             % truth_file
  file2 = arg2;             % diagn_file
  pinfo_out.truth_file = file1;
  pinfo_out.diagn_file = file2;
else
  error('Wrong number of arguments: must be 1 (pinfo) or 2 (file1,file2)')
end

% set this up for later
pinfo_out = setfield(pinfo_out, 'truth_time', [-1,-1]);
pinfo_out = setfield(pinfo_out, 'diagn_time', [-1,-1]);


if ( exist(file1) ~= 2 )
   error(sprintf('(file1) %s does not exist.',file1))
end
if ( exist(file2) ~= 2 )
   error(sprintf('(file2) %s does not exist.',file2))
end

% Get some information from the file1
% dimensions are f1(xxxx), variables are f1{xxxx}
f1 = netcdf(file1);
tmodel      = f1.model(:);
if (isempty(tmodel)) 
   error(sprintf('%s has no ''model'' global attribute.',file1))
end

if ( VarExist(f1,'copy') ) 
   tnum_copies = length(f1('copy')); % determine # of ensemble members
else
   error(sprintf('%s has no ''copy'' dimension.',file1))
end
if ( VarExist(f1,'time') ) 
   tnum_times  = length(f1('time')); % determine # of output times
else
   error(sprintf('%s has no ''time'' dimension.',file1))
end
ttimes = f1{'time'}(:);

[tnum_vars,tdims] = ModelDimension(f1,tmodel);
if ( tnum_vars <= 0 )
   error(sprintf('Unable to determine resolution of %s.',file1))
end
close(f1); 


% Get some information from the file2
f2 = netcdf(file2);
dmodel      = f2.model(:);
if (isempty(dmodel)) 
   error(sprintf('%s has no ''model'' global attribute.',file2))
end
if (VarExist(f2,'copy')) 
   dnum_copies = length(f2('copy')); % determine # of ensemble members
else
   error(sprintf('%s has no ''copy'' dimension.',file2))
end
if (VarExist(f2,'time')) 
   dnum_times  = length(f2('time')); % determine # of output times
else
   error(sprintf('%s has no ''time'' dimension.',file2))
end
dtimes = f2{'time'}(:);
[dnum_vars,ddims] = ModelDimension(f2,dmodel);
if ( dnum_vars <= 0 )
   error(sprintf('Unable to determine resolution of %s.',file2))
end
close(f2); 

% rudimentary bulletproofing
if (strcmp(tmodel,dmodel) ~= 1)
   disp(sprintf('%s has model %s ',file1,tmodel))
   disp(sprintf('%s has model %s ',file2,dmodel))
   error('no No NO ... models must be the same')
end
if (prod(tnum_vars) ~= prod(dnum_vars))
   disp(sprintf('%s has %d state variables',file1,prod(tnum_vars)))
   disp(sprintf('%s has %d state variables',file2,prod(dnum_vars)))
   error('no No NO ... both files must have same shape of state variables.')
end

% if the lengths of the time arrays did not match, this used to be an
% error.  now we call a function to try to find any overlapping regions
% in the time arrays and pass them back up to the called in the pinfo struct.
% they then get used to extract the corresponding hyperslabs of data for
% the matching times.

% construct the pinfo struct in this function
pinfo_out = timearray_intersect(pinfo_out, file1, file2, ttimes, dtimes);

% fail here if the times had nothing in common.
if ( ( pinfo_out.truth_time(1) == -1 ) || ...
     ( pinfo_out.truth_time(2) == -1 ) || ...
     ( pinfo_out.diagn_time(1) == -1 ) || ...
     ( pinfo_out.diagn_time(2) == -1 ))
   disp(sprintf('%s has %d timesteps, from %f to %f', ...
                file1,tnum_times,ttimes(1), ttimes(tnum_times)))  
   disp(sprintf('%s has %d timesteps, from %f to %f', ...
                file2,dnum_times,dtimes(1), dtimes(dnum_times)))  
   error('These files have no timesteps in common')
end


%----------------
% min1,max1 and min2,max2 are the index numbers of the intersection of the
% two input arrays.  -1s in those numbers means no intersection.  1, length()
% means identical (could add a separate flag to simplify the calling code).
function pret = timearray_intersect(pinfo, file1, file2, times1, times2)

% for floating point comparisons, must be within this (single precision)
% roundoff
epsilon = 0.0000001;

% default return; no intersection
pret = pinfo;
pret.truth_file = file1;
pret.diagn_file = file2;
pret.truth_time = [-1,-1];
pret.diagn_time = [-1,-1];

%disp('at start')
%pret

% ensure times are increasing and monotonic, and do they need to be
% a constant delta or not?  compute delta array and validate those match?
% (to within an epsilon with floating pt roundoff)

% check for the no-brainer case - identical time arrays.
% watch out for the floating point compares, and the min/max are probably
% redundant with the (1) and (l) comparisons, but until we put in checks
% for monotonicity, it's a cheap safety check.
l = length(times1);
if (   (length(times1) == length(times2)) ...
    && (abs(min(times1) - min(times2)) < epsilon) ...
    && (abs(max(times1) - max(times2)) < epsilon) ...
    && (times1(1) == times2(1)) ...
    && (times1(l) == times2(l)))
  pret.truth_time = [1,l];
  pret.diagn_time = [1,l];
%disp('equal')
%pret
  return
end

% A is whichever array has the lower min.  this reduces the number of
% cases below we have to check for.
if (min(times1) < min(times2))
  A = times1;
  B = times2;
else
  B = times1;
  A = times2;
end

% precompute the data max, min, lengths using the A,B assignments
% also, if differences are < epsilon, force equality to simplify
% the comparison code below
lenA = length(A);
lenB = length(B);
minA = min(A);
minB = min(B);
maxA = max(A);
maxB = max(B);
if (abs(minA - minB) < epsilon) , minB = minA; , end
if (abs(maxA - minB) < epsilon) , minB = maxA; , end
if (abs(maxA - maxB) < epsilon) , maxB = maxA; , end

% case 1: disjoint regions; simply return here because 
% return struct was initialized to the 'no intersection' case.
if ((minA < minB) && (maxA < minB))  
  return
end

% case 2: B fully contained in A; return corresponding index nums of overlap
% include equal start & end points in this case.
if ((minA <= minB) && (maxB <= maxA)) 
  minI = find(abs(A - minB) < epsilon);
  maxI = find(abs(A - maxB) < epsilon);
  minJ = 1;
  maxJ = lenB;
else
% case 3: partial overlap, minA lower than minB
  minI = find(abs(A - minB) < epsilon);
  maxI = lenA;
  minJ = 1;
  maxJ = find(abs(B - maxA) < epsilon);
end

% now map back to the original input order - this test must match exactly
% the one used initially to assign A and B above.
if (min(times1) < min(times2))
  min1 = minI;
  max1 = maxI;
  min2 = minJ;
  max2 = maxJ;
else
  min1 = minJ;
  max1 = maxJ;
  min2 = minI;
  max2 = maxI;
end

% now put the indices in the return struct and we are done.
pret.truth_time = [min1,max1];
pret.diagn_time = [min2,max2];

% return here


function x = VarExist(ncid,varname)

x = 0;   % false ... assumed not to exist.
variables = var(ncid);
for i=1:length(variables)
   if ( strmatch(name(variables{i}), varname) == 1 )
      x = 1;   % true ... variables exists in the netcdf file.
   end
end


function x = DimExist(ncid,dimname)

x = 0;   % false ... assumed not to exist.
dimensions = dim(ncid);
for i=1:length(dimensions)
   if ( strmatch(name(dimensions{i}), dimname) == 1 )
      x = 1;   % true ... variables exists in the netcdf file.
   end
end


function [x,y] = ModelDimension(ncid,modelname)

x = 0;
y = NaN;

switch lower(modelname)

   case 'cam'
      lonexist = VarExist(ncid,'lon'  );
      latexist = VarExist(ncid,'lat'  );
      lvlexist = VarExist(ncid,'lev');
      if ( latexist && lonexist && lvlexist )
         dnum_lons = prod(size(ncid('lon')));
         dnum_lats = prod(size(ncid('lat')));
         dnum_lvls = prod(size(ncid('lev')));
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];
      end

   case 'fms_bgrid'
      lonexist = VarExist(ncid,'TmpI'  );
      latexist = VarExist(ncid,'TmpJ'  );
      lvlexist = VarExist(ncid,'level');
      if ( latexist && lonexist && lvlexist )
         dnum_lons = prod(size(ncid('TmpI')));
         dnum_lats = prod(size(ncid('TmpJ')));
         dnum_lvls = prod(size(ncid('level')));
         x = 3;
         y = [dnum_lons dnum_lats dnum_lvls];
      end

   case 'simple_advection'
      if ( VarExist(ncid,'loc1d')) 
         y = prod(size(ncid('loc1d')));
	 x = 1;
      end

   otherwise
%     disp(sprintf('working with %s',modelname))
      if ( VarExist(ncid,'StateVariable')) 
         y = prod(size(ncid('StateVariable')));
	 x = 1;
      end

end

