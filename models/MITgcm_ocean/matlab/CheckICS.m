function CheckICS()
% CheckICS 
% 
%

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

%-------------------------------------------------------------------------------
dsize = [256 225 40];
%-------------------------------------------------------------------------------

mitbase = '/ptmp/thoar/MITgcm/ics';
fname   = sprintf('%s/filter_ics',mitbase);
ic00    = rdinit(fname,dsize);
orgmask = find( ic00.data == 0.0 ); % should be dry locations

if ( 1 == 1 )
   
   for i = 1:40
      fname = sprintf('%s/ens_mem_%03d',mitbase,i);
      copy  = rdinit(fname,dsize);
      dry   = copy.data(orgmask);     % should all be dry ... 0.0
      tdiff = sum(dry(:));
      disp(sprintf('copy %02d total is %bx aka %e',i,tdiff,tdiff))
      
      mymask = find(copy.data == 0.0);
      disp(sprintf('mask %i length %06d mask 1 length %06d',i,length(orgmask),length(mymask)))
      if (length(orgmask) ~= length(mymask)) 
         error(sprintf('mask %i length unequal to mask 1 length',i))
      end
   
%     bob = copy.data(:) - ic00.data(:);
%     disp(sprintf('tot diff between ens_mem_%03d and ens_mem_001 %f',i,sum(bob(isfinite(bob)))))

      disp(sprintf('timestamp of ens_mem_%03d %f %f\n',i,copy.time(1),copy.time(2)))
   
   end
   
else


end

%======================================================================

function y = rdinit(fname,dsize)

if (exist(fname,'file'))
  % disp(sprintf('Opening %s',fname))
else
   error('Opening %s',fname)
end

% A DART Initial Conditions file
fid       = fopen(fname,'rb','ieee-be');
h1        = fread(fid,1,'int32');
t         = fread(fid,2,'int32');
h2        = fread(fid,1,'int32');

h3        = fread(fid,    1,'int32');
[x count] = fread(fid,prod(dsize),'float64');
h4        = fread(fid,    1,'int32');
fclose(fid);

if (count ~= prod(dsize))
   error('only read %d of %d items from %s',count,prod(dsize),fname)
end

y.data = reshape(x, dsize);
y.time = t;
