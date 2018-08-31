% check that qice has qualities I think it does qice=enthalpy of ice

%ncf='/glade/scratch/yfzhang/da_test/run/assimilate_ice/cice_prior.r.0010.nc';  % a good restart
ncf='./bad.cice.r.nc'
%ncf='./dart_restart.nc'

aicen=nc_varget(ncf,'aicen');
vicen=nc_varget(ncf,'vicen');
qice=zeros(8,5,384,320);
for lyr=1:8
  qice(lyr,:,:,:)=nc_varget(ncf,['qice00',num2str(lyr)]);
end
aice=squeeze(sum(aicen));
%pcolor(squeeze(aice)); shading flat; colorbar

lyr=4; cat=1;
figure(1); clf
   subplot(211);
   pcolor(squeeze(aicen(cat,:,:)));
   shading flat;
   colorbar

   subplot(212);
   pcolor(squeeze(qice(lyr,cat,:,:)));
   shading flat;
   colorbar

% qice is negative when there is ice or zero when none
% when qice is zero aicen is also zero
% when qice is negative aicen is not zero
% virtually certain qice is not volume weighted, as I thought
i=285:384; cat=3;
x=aicen(cat,i,:); y=qice(lyr,cat,i,:); z=vicen(cat,i,:); xx=aice(i,:);

figure(2); clf
   subplot(221); plot(x(:),y(:),'.'); xlabel('aicen'), ylabel('qice');
   subplot(222); plot(x(:),y(:).*z(:),'.'); xlabel('aicen'), ylabel('qice*vicen');

   subplot(223); plot(xx(:),y(:),'.'); xlabel('aice'), ylabel('qice');
   subplot(224); plot(xx(:),y(:).*z(:),'.'); xlabel('aice'), ylabel('qice*vicen');

figure(3); clf
   k=find(y==0); 
   subplot(211); 
   plot(x(k)); 
   title('aicen where qice is zero, should be zero');

   k=find(y<0);
   subplot(212);
   plot(x(k));
   title('aicen where qice is negative, should not be zero') 

find(x(k)==0), % should be empty

figure(4);clf
   k=find(x==0);
   plot(y(k));
   title('qice where aicen is zero, should be zero')
