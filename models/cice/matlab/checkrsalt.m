% check that sice has qualities I think it does sice=salinity of ice

%ncf='/glade/scratch/yfzhang/da_test/run/assimilate_ice/cice_prior.r.0010.nc';  % a good restart
ncf='/glade/scratch/yfzhang/osse_aggre_hicen/run/osse_aggre_hicen.cice_0019.r.2001-01-02-00000.nc';  % a bad one


%ncf = 'dart_restart.nc'  %after correction
%ncf = 'dart_restart.nc.orig' %before correction

aicen=nc_varget(ncf,'aicen');
vicen=nc_varget(ncf,'vicen');
sice=zeros(8,5,384,320);
for lyr=1:8
  sice(lyr,:,:,:)=nc_varget(ncf,['sice00',num2str(lyr)]);
end
aice=squeeze(sum(aicen));
%pcolor(squeeze(aice)); shading flat; colorbar

lyr=4; cat=1;
figure(1); clf;
subplot(221);pcolor(squeeze(aicen(cat,:,:))); shading flat; colorbar
subplot(223);pcolor(squeeze(sice(lyr,cat,:,:))); shading flat; colorbar

subplot(222); plot(squeeze(sice(:,4,360,180)),1:8); % fresher at top
set(gca,'ydir','rev'); ylabel('layer'); xlabel('salinity');
title('older/thicker ice, fresher at top')
subplot(224); plot(squeeze(sice(:,1,360,180)),1:8); 
title('younger/thinner ice, saltier at top & bottom')
set(gca,'ydir','rev'); ylabel('layer'); xlabel('salinity');

% sice is positive when there is ice or zero when none
% when sice is zero aicen is also zero
% when sice is positive aicen is not zero
% virtually certain sice is not volume weighted
i=285:384; cat=1;
x=aicen(cat,i,:); y=sice(lyr,cat,i,:); z=vicen(cat,i,:); xx=aice(i,:);
figure(2); clf; subplot(221); plot(x(:),y(:),'.'); xlabel('aicen'), ylabel('sice');
subplot(222); plot(x(:),y(:).*z(:),'.'); xlabel('aicen'), ylabel('sice*vicen');

subplot(223); plot(xx(:),y(:),'.'); xlabel('aice'), ylabel('sice');
subplot(224); plot(xx(:),y(:).*z(:),'.'); xlabel('aice'), ylabel('sice*vicen');

figure(3); clf
k=find(y==0); subplot(211); plot(x(k)); title('aicen where sice is zero, should be zero');
k=find(y>0); subplot(212);  plot(x(k)); title('aicen where sice is positive, should not be zero') 
find(x(k)==0), % should be empty

figure(4); clf
k=find(x==0);plot(y(k));title('sice where aicen is zero, should be zero');
