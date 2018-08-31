% check that qsno has qualities I think it does qsno=enthalpy of snow

%ncf='/glade/scratch/yfzhang/da_test/run/assimilate_ice/cice_prior.r.0010.nc';  % a good restart
ncf='/glade/scratch/yfzhang/osse_aggre_hicen/run/osse_aggre_hicen.cice_0019.r.2001-01-02-00000.nc';  % a bad one


%ncf = './dart_restart.nc'          %after correction
%ncf = './dart_restart.nc.orig'     %before correction

aicen=nc_varget(ncf,'aicen');
vicen=nc_varget(ncf,'vicen');
qsno=zeros(3,5,384,320);
for lyr=1:3
  qsno(lyr,:,:,:)=nc_varget(ncf,['qsno00',num2str(lyr)]);
end
aice=squeeze(sum(aicen));
%pcolor(squeeze(aice)); shading flat; colorbar

lyr=2; cat=1;
figure(1); clf
subplot(211);pcolor(squeeze(aicen(cat,:,:))); shading flat; colorbar
subplot(212);pcolor(squeeze(qsno(lyr,cat,:,:))); shading flat; colorbar


% qsno is negative when there is ice or zero when none
% when qsno is zero aicen is also zero
% when qsno is negative aicen is not zero
% virtually certain qsno is not volume weighted, as I thought
i=285:384; cat=1;
x=aicen(cat,i,:); y=qsno(lyr,cat,i,:); z=vicen(cat,i,:); xx=aice(i,:);
figure(2); clf
subplot(221); plot(x(:),y(:),'.'); xlabel('aicen'), ylabel('qsno');
subplot(222); plot(x(:),y(:).*z(:),'.'); xlabel('aicen'), ylabel('qsno*vicen');

subplot(223); plot(xx(:),y(:),'.'); xlabel('aice'), ylabel('qsno');
subplot(224); plot(xx(:),y(:).*z(:),'.'); xlabel('aice'), ylabel('qsno*vicen');

figure(3); clf
k=find(y==0); subplot(211); plot(x(k)); title('aicen where qsno is zero, should be zero');
k=find(y<0); subplot(212);  plot(x(k)); title('aicen where qsno is negative, should not be zero') 
find(x(k)==0), % should be empty

figure(4); clf
k=find(x==0);  plot(y(k)); title('qsno where aicen is zero, should be zero');

