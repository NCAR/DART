
ncf='cice_out.r.0001.nc';
%ncf='cice_in.r.0002.nc';
aicenout=nc_varget(ncf,'aicen');
aout=squeeze(sum(aicenout));
ncf='cice_in.r.0001.nc';
aicenin=nc_varget(ncf,'aicen');
ain=squeeze(sum(aicenin));
diff=ain-aout;

subplot(211); pcolor(ain+diff*.3e8); shading flat; axis([120 250 280 380])
subplot(212); pcolor(ain-aout); shading flat; axis([120 250 280 380])
return

subplot(211); pcolor(squeeze(aicenout(1,:,:))); shading flat;
subplot(212); pcolor(squeeze(aicenout(1,:,:)-aicenin(1,:,:))); shading flat;

return

ncf='cice.r.nc';

aicenin=nc_varget(ncf,'aicen');

aice=squeeze(sum(aicenin));
pcolor(squeeze(aice)); shading flat; colorbar

