function fit_ens_spread_time
%
p=load('Tanl_times_level.dat');
level=p(1);

   figure(1); clf
%  orient landscape
   orient portrait
   wysiwyg

   ges   = sprintf('Tges_times_%04dmb.dat',level);
   anl   = sprintf('Tanl_times_%04dmb.dat',level);

   region = 'NH';
   top   = sprintf('T ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,1)
   
   region = 'SH';
   top   = sprintf('T ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,2)
   
   region = 'TR';
   top   = sprintf('T ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,3)
   
   region = 'NA';
   top   = sprintf('T ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,4)
   
   % Now for the Winds
   figure(2); clf
%  orient landscape
   orient portrait
   wysiwyg
   
   ges   = sprintf('Wges_times_%04dmb.dat',level);
   anl   = sprintf('Wanl_times_%04dmb.dat',level);

   region = 'NH';
   top   = sprintf('Wind ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,1)
   
   region = 'SH';
   top   = sprintf('Wind ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,2)
   
   region = 'TR';
   top   = sprintf('Wind ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,3)
   
   region = 'NA';
   top   = sprintf('Wind ens spread %04dmb %s',level,region);
   Myplot(ges,anl,top,4)

   str = '-dpsc ';

   print(1,str,'t_ens_spread_time.ps');
   print(2,str,'w_ens_spread_time.ps');

% end



function Myplot(file1,file2,top,region)

p=load(file1);
a=load(file2);
xp=p(:,1);
xa=a(:,1);
count=3+(region-1)*3;
yp_spread=p(:,count);
ya_spread=a(:,count);
%
subplot(2,2,region)
plot(xp,yp_spread,'c+-',xa,ya_spread,'ro:')
grid
%if region > 2 
xlabel('Time interval', 'fontsize', 10) ;
%end
ylabel('RMSE ', 'fontsize', 10)
%title({'Fit to RAOBS, ensemble spread',file1,file2}, ...
%      'fontsize', 10,'Interpreter','none')
title({top}, ...
      'fontsize', 10,'Interpreter','none')
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
