%module load matlab
%matlab -nodisplay -nodesktop -nosplash 

% DART $Id$
% CREDIT: Alexey Morozov

clear
clc
close all
format compact

figure(1); clf
ne=3;
nx=3; %inf
ny=1; %aw

type='mult';
ispost='T';

for k=1:ny
  for j=1:nx
  cd(['/nobackup/morozova/' type '_' ispost '_' num2str(k) '_' num2str(j) '/gitm/work/']) 
    subplot(ny,nx,j+nx*(k-1))
    hold on
    for i=1:ne 
  r=0;
try
  r=load([ num2str(i) '.dat']);
stairs(r,'color',rand(1,3))
%  xlim([1 1440/30])
end
  end
hold off
end
end

cd ~/

legend('1','2','3')

%% printing
set(gcf,'papertype','usletter')
  set(gcf,'paperorientation','landscape')
  set(gcf,'paperposition',[.25 .25 10.5 8])
  print(gcf, '-dpdf',[ type '_' ispost ])

  unix(['scp ~/' type '_' ispost '.pdf login.engin.umich.edu:Public/html/um']);

disp(['umich.edu/~morozova/um/' type '_' ispost '.pdf']  )

cd ~/gitm11/matlab

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
