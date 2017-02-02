function VerifyNetCDF(varname)
%% VerifyNetCDF compares the fortran direct-access input files with the variables
% in the True_State.nc file. This script was used to ensure the state variable was
% being read properly and that we could parse it correctly ... and that the netCDF
% routines were working as expected.
% 
% VerifyNetCDF('T')

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

varname = upper(varname);
mitbase = '/fs/image/home/nancy/subversion/trunk/models/MITgcm_ocean/data2/';

switch varname
case {'S','T','U','V'}
   mitO  = rdmds(sprintf('%s/%s.0000040992',mitbase,varname));
   [nx ny nz] = size(mitO);
otherwise
   mitO  = rdmds(sprintf('%s/Eta.0000040992',mitbase));
   [nx ny] = size(mitO);
   nz = 1;
   varname = 'Eta';
end

levels = getnc('True_State.nc','ZC');

% uncomfortable assumption about missing value flag
% The DART routines make the same uncomfortable assumption.

inds   = find(mitO == 0.0);
mitO(inds) = NaN;
   
%  fid = fopen('perfect_ics','rt');
%  [validtime, count] = fscanf(fid,'%d',2);
%  x     = fscanf(fid,'%f');
%  fclose(fid)
%  
%  mysize = prod(size(s)) + prod(size(t)) + prod(size(u)) + prod(size(v)) + prod(size(ssh));
%  psize = length(x);
%  
%  disp(sprintf('MIT size is %d, DART size is %d',mysize,psize))

%  s_end = prod(size(s));
%  ICS = reshape(x(1:s_end),size(s));

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

for i = 1:nz

   figure(1); clf; colormap(gauss3);

   switch varname
   case {'S','T','U','V'}
      mit = squeeze(mitO(:,:,i));
      corner = [-1 -1 i -1 -1];
      endpnt = [-1 -1 i -1 -1];
   otherwise
      mit = mitO(:,:);
      corner = [-1 i -1 -1];
      endpnt = [-1 i -1 -1];
   end

   dart = getnc('True_State.nc',varname,corner,endpnt,-1,-2);

   datmat = (dart - mit);

   inds = isfinite(mit);
   mmin = min(mit(inds));
   mmax = max(mit(inds));
   nmitzeros = length(inds);
   b = datmat(inds);

   inds = isfinite(dart);
   dmin = min(dart(inds));
   dmax = max(dart(inds));
   ndartzeros = length(inds);

   if (ndartzeros ~= nmitzeros) 
      disp(sprintf('WARNING: mit has %d zeros dart has %d zeros',nmitzeros,ndartzeros))
   end
   if (mmin ~= dmin) 
      disp(sprintf('WARNING: mit minimum %f dart minimum %d',mmin,dmin))
   end
   if (mmax ~= dmax) 
      disp(sprintf('WARNING: mit maximum %f dart maximum %d',mmax,dmax))
   end

   subplot(2,2,1)
   imagesc(mit'); set(gca,'YDir','normal')
   title({sprintf('Original : level %d of %d   %f',i,nz,levels(i)),...
          sprintf('min %f max %f',mmin,mmax)})
   colorbar

   subplot(2,2,2)
   imagesc(dart'); set(gca,'YDir','normal')
   title(sprintf('netcdf min %f max %f',dmin,dmax))
   colorbar

   subplot(2,2,3)
   imagesc(datmat'); set(gca,'YDir','normal')
   title({'difference (netcdf - original)', ...
         sprintf('min %f max %f',min(b), max(b))})
   colorbar

   subplot(2,2,4)
   hist(b,50);
   title({'difference (netcdf - original)', ...
         sprintf('min %f max %f',min(b), max(b))})

   % disp('Pausing - hit any key to continue')
   disp('Pausing for one second')
   pause(1.0)

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
