function fit_ens_mean_time

p=load('Tanl_times_level.dat');
level=p(1);

   figure(1); clf
%  orient landscape
   orient portrait
   wysiwyg

   ges   = sprintf('Tges_times_%04dmb.dat',level);
   region = 'NH';
   top   = sprintf('T obs number            %04dmb %s',level,region);
   Myplot(ges,top,1)
   
   region = 'SH';
   top   = sprintf('T obs number            %04dmb %s',level,region);
   Myplot(ges,top,2)
   
   region = 'TR';
   top   = sprintf('T obs number            %04dmb %s',level,region);
   Myplot(ges,top,3)
   
   region = 'NA';
   top   = sprintf('T obs number            %04dmb %s',level,region);
   Myplot(ges,top,4)
   
   % Now for the Winds
   figure(2); clf
%  orient landscape
   orient portrait
   wysiwyg
   
   ges   = sprintf('Wges_times_%04dmb.dat',level);
   region = 'NH';
   top   = sprintf('Wind obs number            %04dmb %s',level,region);
   Myplot(ges,top,1)
   
   region = 'SH';
   top   = sprintf('Wind obs number            %04dmb %s',level,region);
   Myplot(ges,top,2)
   
   region = 'TR';
   top   = sprintf('Wind obs number            %04dmb %s',level,region);
   Myplot(ges,top,3)
   
   region = 'NA';
   top   = sprintf('Wind obs number            %04dmb %s',level,region);
   Myplot(ges,top,4)

   str = '-dpsc ';

   print(1,str,'t_obs_num_time.ps');
   print(2,str,'w_obs_num_time.ps');

% end

function Myplot(file1,top,region)

p=load(file1);
xp=p(:,1);
count=4+(region-1)*3;
%count
yp_num=p(:,count);
%
subplot(2,2,region)
plot(xp,yp_num,'c+-')
grid
%if region > 2 
xlabel('Time interval', 'fontsize', 10) ;
%end
%ylabel('RMSE ', 'fontsize', 10)
%title({'Fit to RAOBS, ensemble mean',file1,file2}, ...
%      'fontsize', 10,'Interpreter','none')
title({top}, ...
      'fontsize', 10,'Interpreter','none')
%legend('guess', 'analysis')
%h=legend;
%legend(h,'boxoff')
