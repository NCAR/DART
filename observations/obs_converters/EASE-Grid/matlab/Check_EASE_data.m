%
%
%

%%-------------------------------------------------------------------------------
% Read what I wrote from the fortran program
%--------------------------------------------------------------------------------
fid = fopen('/glade/p/work/thoar/DART/Tb/observations/AMSR-E/work/sanity_check.ieee','rb');

bob = fread(fid,1,'int32');
nrows = fread(fid,1,'int32') 
ncols = fread(fid,1,'int32') 
bob = fread(fid,1,'int32');

bob = fread(fid,1,'int32');
datmat = fread(fid,[nrows,ncols],'int32');
bob = fread(fid,1,'int32');
fclose(fid);

%%-------------------------------------------------------------------------------
% Read the original data
%--------------------------------------------------------------------------------
AMSR = '/glade/p/image/Observations/land/EASE_Grid/2004_north/ID2r3-AMSRE-NL2004004A.v03.89V';
fid  = fopen(AMSR,'rb');
Tb   = fread(fid,[nrows,ncols],'uint16');
fclose(fid)

%%-------------------------------------------------------------------------------
% Compare the two
%--------------------------------------------------------------------------------

figure(1); clf
imagesc(datmat);
set(gca,'Ydir','normal');
colorbar
title('Fortran')


figure(2); clf
imagesc(Tb);
set(gca,'Ydir','normal');
colorbar
title('AMSR-E')

mydiff = datmat - Tb;
[min(    Tb(:)) max(    Tb(:))] 
[min(datmat(:)) max(datmat(:))] 
[min(mydiff(:)) max(mydiff(:))] 

if ( 1 == 2 )
   figure(3); clf
   fname         = '/Users/thoar/svn/DART/Tb/observations/AMSR-E/work/obs_epoch_001.nc';
   ObsTypeString = 'AMSRE_BRIGHTNESS_T';
   region        = [0 360 0 90 -Inf Inf];
   CopyString    = 'observation';
   QCString      = 'Data QC';
   maxgoodQC     = 2;
   verbose       = 1;   % anything > 0 == 'true'
   twoup         = 0;   % anything > 0 == 'true'
    
   bob = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, maxgoodQC, verbose, twoup);
end

