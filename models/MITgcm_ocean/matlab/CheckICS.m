function CheckICS()
%% CheckICS 
% 
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%-------------------------------------------------------------------------------
dsize   = [256 225 40];
ncopies = 40;
%-------------------------------------------------------------------------------

if ( 1 == 2 )
   
   mitbase = '/ptmp/thoar/MITgcm/ics';
   fname   = sprintf('%s/filter_ics',mitbase);
   ic00    = rdinit(fname,dsize);
   orgmask = find( ic00.data == 0.0 ); % should be dry locations

   for i = 1:ncopies
      fname = sprintf('%s/ens_mem_%03d',mitbase,i);
      copy  = rdinit(fname,dsize);
      dry   = copy.data(orgmask);     % should all be dry ... 0.0
      tdiff = sum(dry(:));
      disp(sprintf('copy %02d total is %bx aka %e',i,tdiff,tdiff))
      
      mymask = find(copy.data == 0.0);
      disp(sprintf('mask %i length %06d mask 1 length %06d',i,length(orgmask),length(mymask)))
      if (length(orgmask) ~= length(mymask)) 
         error('mask %i length unequal to mask 1 length',i)
      end
   
%     bob = copy.data(:) - ic00.data(:);
%     disp(sprintf('tot diff between ens_mem_%03d and ens_mem_001 %f',i,sum(bob(isfinite(bob)))))

      disp(sprintf('timestamp of ens_mem_%03d %f %f\n',i,copy.time(1),copy.time(2)))
   
   end
   
else

   mitbase = '/ptmp/thoar/MITgcm/ics';
   fname   = sprintf('%s/filter_ics',mitbase);
   ncopies = 40;

   mitbase = '/ptmp/thoar/MITgcm/obs_gliders';
   fname   = sprintf('%s/filter_restart',mitbase);
   ncopies = 20;

   fid     = fopen(fname,'rb','ieee-be');

   for i = 1:ncopies

      % Read the time stamp for each ensemble member
      h1        = fread(fid,1,'int32');
      T         = fread(fid,2,'int32');
      h2        = fread(fid,1,'int32');
      
      if ( h1 ~= h2 )
          error('copy %02d time record lengths %d = %d',i,h1,h2)
      else
          disp(sprintf('Copy %02d timestamp is days,seconds %d,%d',i,T(2),T(1)))
      end

      % Parse the state vector
      h3  = fread(fid,    1,'int32');
      s   = fread(fid,prod(dsize),'float64');
      t   = fread(fid,prod(dsize),'float64');
      u   = fread(fid,prod(dsize),'float64');
      v   = fread(fid,prod(dsize),'float64');
      eta = fread(fid,prod(dsize(1:2)),'float64');
      h4  = fread(fid,    1,'int32');
      if ( h3 ~= h4 )
          error('copy %02d data record lengths %d = %d',i,h3,h4)
      end

      % Plot salinity for each ensemble member
      Y         = reshape(s, dsize);
      Y(Y==0.0) = NaN;

      for ilev = 1:size(Y,3)
         imagesc(squeeze(Y(:,:,ilev)'))
         set(gca,'YDir','normal')
         title(sprintf('S copy %02d level %03d day %d sec %d\n',i,ilev,T(2),T(1)))
         colorbar;
         display('Pausing, hit any key to continue ...')
         pause
      end

      % Plot temperature for each ensemble member
      Y         = reshape(t, dsize);
      Y(Y==0.0) = NaN;

      for ilev = 1:size(Y,3)
         imagesc(squeeze(Y(:,:,ilev)'))
         set(gca,'YDir','normal')
         title(sprintf('T copy %02d level %03d day %d sec %d\n',i,ilev,T(2),T(1)))
         colorbar;
         display('Pausing, hit any key to continue ...')
         pause
      end

      % Read the U current component for each ensemble member
      Y         = reshape(u, dsize);
      Y(Y==0.0) = NaN;

      for ilev = 1:size(Y,3)
         imagesc(squeeze(Y(:,:,ilev)'))
         set(gca,'YDir','normal')
         title(sprintf('U copy %02d level %03d day %d sec %d\n',i,ilev,T(2),T(1)))
         colorbar;
         display('Pausing, hit any key to continue ...')
         pause
      end

      % Read the V current component for each ensemble member
      Y         = reshape(v, dsize);
      Y(Y==0.0) = NaN;

      for ilev = 1:size(Y,3)
         imagesc(squeeze(Y(:,:,ilev)'))
         set(gca,'YDir','normal')
         title(sprintf('V copy %02d level %03d day %d sec %d\n',i,ilev,T(2),T(1)))
         colorbar;
         display('Pausing, hit any key to continue ...')
         pause
      end

      % Read the sea surface height component for each ensemble member
      Y         = reshape(eta, dsize(1:2));
      Y(Y==0.0) = NaN;

      imagesc(Y')
      set(gca,'YDir','normal')
      title(sprintf('SSH copy %02d level %03d day %d sec %d\n',i,ilev,T(2),T(1)))
      colorbar;
      display('Pausing, hit any key to continue ...')
      pause

   end

   fclose(fid);
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
if ( h1 ~= h2 )
    error('time record lengths %d = %d',h1,h2)
end

h3        = fread(fid,    1,'int32');
[x count] = fread(fid,prod(dsize),'float64');
h4        = fread(fid,    1,'int32');
if ( h3 ~= h4 )
    error('data record lengths %d = %d',h3,h4)
end
fclose(fid);

if (count ~= prod(dsize))
   error('only read %d of %d items from %s',count,prod(dsize),fname)
end

y.data = reshape(x, dsize);
y.time = t;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
