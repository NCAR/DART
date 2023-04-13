function [] = plot_diagnostic(yearstr,Faxa_swndr,Faxa_swndf,Faxa_swvdr,Faxa_swvdf, ...
                              Faxa_rainl,Faxa_rainc,Faxa_snowl,Faxa_snowc,Sa_tbot, ...
                              Sa_shum,Sa_u,Sa_v,Sa_pbot,Faxa_lwdn, ...
                              Faxa_swndr_adjust,Faxa_swndf_adjust,Faxa_swvdr_adjust ...
                              Faxa_swvdf_adjust,Faxa_rainl_adjust,Faxa_snowl_adjust ...
                              Sa_tbot_adjust,Sa_shum_adjust,Sa_u_adjust,Sa_v_adjust ...
                              Sa_pbot_adjust,Faxa_lwdn_adjust, ...
                              sw_1hr,ppt_1hr,ta_1hr,q_1hr,wind_1hr,ps_1hr,lw_1hr ...
                              ppt_3hr,ta_3hr);
% function plot_diagnostic: Full diagnostic package for CAM bias correction
% Output: empty

     % Comparing original CAM6 and Tower met forcing
     % Check for time synchronization, ensemble adjustment, snow/rain partitioning etc. 

       figure(1)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       j=plot(mean(Faxa_swndr,1)+ ...
              mean(Faxa_swndf,1)+ ...
              mean(Faxa_swvdr,1)+ ...
              mean(Faxa_swvdf,1), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(sw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);

       axes(ha(2));
       j=plot(repelem(mean(Faxa_rainl,1),1,3)'+ ...
              repelem(mean(Faxa_rainc,1),1,3)'+ ...
              repelem(mean(Faxa_snowl,1),1,3)'+ ...
              repelem(mean(Faxa_snowc,1),1,3)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)

       axes(ha(3));
       j=plot(repelem(mean(Sa_tbot,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(4));
       j=plot(repelem(mean(Sa_shum,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(q_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Specific Humidity (kg kg^{-1})','FontSize', fontpt)
       linkaxes(ha(1:4), 'x');




       figure(2)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(3,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       j=plot(mean(Sa_u,1), ...
              '--', 'color',rgb('darksalmon'), 'linewidth', mean_width); hold on
       k=plot(mean(Sa_v,1), ...
              '-', 'color',rgb('darksalmon'), 'linewidth', mean_width);
       l=plot( (mean(Sa_u,1).^2+ ...
               mean(Sa_v,1).^2).^0.5, ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width);

       m=plot(wind_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k,l,m], 'Zonal, CAM6 Ens Mean', 'Meridional, CAM6 Ens Mean', ...
                     'Total, CAM6 Ens Mean','Tower Met','FontSize',fontpt);

       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);


       axes(ha(2));
       j=plot(repelem(mean(Sa_tbot,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(3));
       j=plot(repelem(mean(Faxa_rainl,1),1,3)'+ ...
              repelem(mean(Faxa_rainc,1),1,3)'+ ...
              repelem(mean(Faxa_snowl,1),1,3)'+ ...
              repelem(mean(Faxa_snowc,1),1,3)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)

       linkaxes(ha(1:3), 'x');


       figure(3)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(2,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));
       j=plot(repelem(mean(Sa_pbot,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ps_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Pressure (Pa)','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);

       axes(ha(2));
       j=plot(repelem(mean(Faxa_lwdn,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(lw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('LW radiation (W m^{-1})','FontSize', fontpt)
       linkaxes(ha(1:2), 'x');

display('  ')
display('Finished tower CAM6 diagnostic plots. Press enter to proceed to adjusted CAM6 plots..')
display('  ')
pause
close all


     % For plotting purposes only, make these 3 hour products into 1 hour products for easier comparison 
     % Matrices (80, 2920) --> (80, 8760); 

       figure(4)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3)/2 scrsz(4)]);
       fontpt=16;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       ax1=plot(Faxa_swndr(:,188*24:195*24)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust(:,188*24:195*24)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       j=plot(mean(Faxa_swndr(:,188*24:195*24),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Faxa_swndr_adjust(:,188*24:195*24),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt);
       ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
       title(['CAM6 scaling against tower NR1 Met: July ', yearstr] , 'FontSize',fontpt);
       xticks([24:24:24*18]);
       set(ha(1),'xticklabel',[])

       axes(ha(2));

       ax2=plot(repelem(Faxa_rain(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Faxa_rainl_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Faxa_rain(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Faxa_rainl_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Rain (mm  s^{-1})','FontSize', fontpt)
       xticks([24:24:24*18]);
       set(ha(2),'xticklabel',[])

       axes(ha(3));

       ax3=plot(repelem(Sa_tbot(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_tbot_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_tbot(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_tbot_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Temperature (K)','FontSize', fontpt)
       xticks([24:24:24*18]);
       set(ha(3),'xticklabel',[])

       axes(ha(4));

       ax4=plot(repelem(Sa_shum(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_shum_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_shum(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_shum_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
       xticks([24:24:24*18]);
       linkaxes(ha(1:4), 'x');


       figure(5)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       plot(Faxa_swndr(:,32*24:45*24)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust(:,32*24:45*24)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       j=plot(mean(Faxa_swndr(:,32*24:45*24),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Faxa_swndr_adjust(:,32*24:45*24),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt);
       ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
       title(['CAM6 scaling against tower NR1 Met: February ', yearstr] , 'FontSize',fontpt);
       xticks([24:24:24*13]);

       axes(ha(2));

       plot(repelem(Faxa_snow(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Faxa_snowl_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Faxa_snow(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Faxa_snowl_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Snow (mm  s^{-1})','FontSize', fontpt)
       xticks([24:24:24*13]);

       axes(ha(3));

       plot(repelem(Sa_tbot(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_tbot_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_tbot(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_tbot_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Temperature (K)','FontSize', fontpt)
       xticks([24:24:24*13]);
       axes(ha(4));

       plot(repelem(Sa_shum(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_shum_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_shum(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_shum_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
       xticks([24:24:24*13]);

       linkaxes(ha(1:4), 'x');


       figure(6)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(5,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));

       j=plot(mean(Faxa_swndr,1)+ ...
              mean(Faxa_swndf,1)+ ...
              mean(Faxa_swvdr,1)+ ...
              mean(Faxa_swvdf,1), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(sw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       title(['CAM6 scaling against tower NR1 Met ', yearstr] , 'FontSize',fontpt);
       ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)


       axes(ha(2));

       ax1=plot(Faxa_swndr', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       j=plot(mean(Faxa_swndr,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Faxa_swndr_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt);
       ylabel('Near-IR, Direct (W m^{-2})','FontSize', fontpt)

       axes(ha(3));

       ax2=plot(Faxa_swndf', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndf_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_swndf,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_swndf_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Near-IR, Diffuse (W m^{-2})','FontSize', fontpt)

       axes(ha(4));

       ax3=plot(Faxa_swvdr', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swvdr_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_swvdr,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_swvdr_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Visible, Direct (W m^{-2})','FontSize', fontpt)

       axes(ha(5));

       ax4=plot(Faxa_swvdf', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swvdf_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_swvdf,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_swvdf_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Visible, Diffuse (W m^{-2})','FontSize', fontpt)
       linkaxes(ha(1:5), 'x');


       figure(7)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));
       j=plot(mean(Faxa_rainl,1)'+ ...
              mean(Faxa_rainc,1)'+ ...
              mean(Faxa_snowl,1)'+ ...
              mean(Faxa_snowc,1)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_3hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       title(['CAM6 scaling against tower NR1 Met ', yearstr] , 'FontSize',fontpt);
       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)
       ylim([0 2*10^-3])

       axes(ha(2));
       j=plot(mean(Sa_tbot,1), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_3hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(3));
       plot(Faxa_rain', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_rainl_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_rain,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_rainl_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Rain (mm  s^{-1})','FontSize', fontpt)
       ylim([0 2*10^-3])

       axes(ha(4));

       plot(Faxa_snow', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_snowl_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_snow,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_snowl_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Snow (mm  s^{-1})','FontSize', fontpt)
       ylim([0 2*10^-3])
       linkaxes(ha(1:4), 'x');

       figure(8)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(5,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));

       j=plot(mean(Sa_u,1), ...
              '--', 'color',rgb('darksalmon'), 'linewidth', mean_width); hold on
       k=plot(mean(Sa_v,1), ...
              '-', 'color',rgb('darksalmon'), 'linewidth', mean_width);
       l=plot( (mean(Sa_u,1).^2+ ...
               mean(Sa_v,1).^2).^0.5, ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width);

       m=plot(wind_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k,l,m], 'Zonal, CAM6 Ens Mean', 'Meridional, CAM6 Ens Mean', ...
                     'Total, CAM6 Ens Mean','Tower Met','FontSize',fontpt);

       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);
       axes(ha(2));

       ax1=plot(Sa_u', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_u_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);

       j=plot(mean(Sa_u,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Sa_u_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt);
       ylabel('Zonal wind (m s^{-1})','FontSize', fontpt)
       axes(ha(3));

       ax2=plot(Sa_v', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_v_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Sa_v,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Sa_v_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Meridional wind (m s^{-1})','FontSize', fontpt)
       axes(ha(4));

       ax3=plot(repelem(Sa_pbot,1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_pbot_adjust,1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_pbot,1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_pbot_adjust,1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Pressure (Pa)','FontSize', fontpt)

       axes(ha(5));

       ax4=plot(repelem(Faxa_lwdn,1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Faxa_lwdn_adjust,1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Faxa_lwdn,1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Faxa_lwdn_adjust,1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Longwave Radiation (W m^{-2})','FontSize', fontpt)
       linkaxes(ha(1:5), 'x');

display('   ')
display(['Adjusted CAM6 plots and diagnostics completed.  Press enter to write ' yearstr ' CAM6 adjusted files to netcdf'])
display('   ')
pause
close all
	

end
