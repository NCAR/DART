function [oplus, cg_height, zkm] = binary_vs_netcdf(binaryfile,file2)
%% function to explore the discontinuity in oplus.
% longitude index 16 has a discontinuity along the latitude dimension for 
% the northern hemisphere high latitudes that is not present at the other
% longitudes. There is a problem with the heights ... the last 4 heights
% are perfect zeros.
%    
%    binaryfile = '../data/dart.oplus2.bin';
%    file2 = '../data/DATA.ionos2.nc';
%    [oplus, height] = binary_vs_netcdf(binaryfile, file2);
%

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% open(33,file='dart.oplus2.bin',form='unformatted')
% write(33)size(d13d,1),size(d13d,2),size(d13d,3)
% write(33)d13d(:,:,:,1)
% write(33)size(zkm)
% write(33)zkm
% close(33)

ggcm = fopen(binaryfile,'r');
byte = fread(ggcm,1,'int32');
dim1 = fread(ggcm,1,'int32');
dim2 = fread(ggcm,1,'int32');
dim3 = fread(ggcm,1,'int32');
byte = fread(ggcm,1,'int32');

byte = fread(ggcm,1,'int32');
gdat = fread(ggcm,dim1*dim2*dim3,'float64');
byte = fread(ggcm,1,'int32');

byte = fread(ggcm,1,'int32');
dim1 = fread(ggcm,1,'int32');
byte = fread(ggcm,1,'int32');

byte = fread(ggcm,1,'int32');
zkm  = fread(ggcm,dim1,'float64');
byte = fread(ggcm,1,'int32');
stat = fclose(ggcm);

datmat = reshape(gdat,dim1,dim2,dim3);

ncid      = netcdf.open(file2,'NC_NOWRITE');
cg_lat    = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'cg_lat'));
cg_lon    = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'cg_lon'));
cg_height = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'cg_height'));

fprintf('cg_lon(15) is %s\n',cg_lon(15))
fprintf('cg_lon(16) is %s\n',cg_lon(16))
fprintf('cg_lon(17) is %s\n',cg_lon(17))

varid = netcdf.inqVarID(ncid,'oplus');
oplus = netcdf.getVar(ncid,varid);
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,varid);

for idim=1:length(dimids)
    [dimname, dimlen] = netcdf.inqDim(ncid,dimids(idim));
    dimstring{idim} = sprintf('%s has length %d',dimname,dimlen);
end
netcdf.close(ncid)

size(datmat)
size(oplus)

for ilon = 1:dim3

   fprintf('longitude %d\n',ilon)

   figure(1); orient tall;
   subplot(2,1,1)
   imagesc(squeeze(datmat(:,:,ilon)));
   
   title({sprintf('%s longitude %d',binaryfile, ilon), sprintf('%d -x- %d -x- %d',dim1,dim2,dim3)})

   ylabel(sprintf('dimension %d (length=%d)',1,dim1))
   xlabel(sprintf('dimension %d (length=%d)',2,dim2))
   axis image
   wysiwyg
   set(gca,'FontSize',20)
   colorbar('eastoutside')
   
   subplot(2,1,2)
   imagesc(squeeze(oplus(:,:,ilon)));

   axis image
   set(gca,'FontSize',20)
   colorbar('eastoutside')
   title('the netCDF equivalent')
   xlabel('index')
   ylabel('index')
   
%   subplot(3,1,3)
%   imagesc(cg_lat,cg_height,squeeze(oplus(:,:,ilon)));
%   axis image
%   set(gca,'Ydir','normal')
%   % set(gca,'FontSize',20)
%   colorbar('eastoutside')
%   title('the netCDF equivalent')
%   xlabel(dimstring{2},'Interpreter','none')
%   ylabel(dimstring{1},'Interpreter','none') 

    fname = sprintf('oplus_lonindex_%d',ilon);
    print(gcf,'-dpdf',fname)
 
   disp('Pausing, hit any key to continue ...')
   pause(0.1)

end



% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

