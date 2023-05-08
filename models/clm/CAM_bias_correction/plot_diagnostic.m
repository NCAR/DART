function [] = plot_diagnostic(yearstr,Faxa_swndr,Faxa_swndf,Faxa_swvdr,Faxa_swvdf, ...
    Faxa_rainl,Faxa_rainc,Faxa_snowl,Faxa_snowc,Sa_tbot, ...
    Sa_shum,Sa_u,Sa_v,Sa_pbot,Faxa_lwdn,Faxa_rain,Faxa_snow, ...
    Faxa_swndr_adjust,Faxa_swndf_adjust,Faxa_swvdr_adjust, ...
    Faxa_swvdf_adjust,Faxa_rainl_adjust,Faxa_snowl_adjust, ...
    Sa_tbot_adjust,Sa_shum_adjust,Sa_u_adjust,Sa_v_adjust, ...
    Sa_pbot_adjust,Faxa_lwdn_adjust,n, ...
    sw_1hr,ppt_1hr,ta_1hr,q_1hr,wind_1hr,ps_1hr,lw_1hr, ...
    ppt_3hr,ta_3hr,wind_6hr,ta_6hr,ppt_6hr,sw_6hr);
% function plot_diagnostic: Full diagnostic package for CAM bias correction
% Input : CAM, corrected CAM, and tower met. Final 4 variables are optional (CAM4 only)
% Output: empty

% Comparing original CAM and Tower met forcing
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
switch n
    case 'CAM6'
        j=plot(mean(Faxa_swndr,1)+ ...
            mean(Faxa_swndf,1)+ ...
            mean(Faxa_swvdr,1)+ ...
            mean(Faxa_swvdf,1), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
    case 'CAM4'
        j=plot(repelem(mean(Faxa_swndr,1),1,6)'+ ...
            repelem(mean(Faxa_swndf,1),1,6)'+ ...
            repelem(mean(Faxa_swvdr,1),1,6)'+ ...
            repelem(mean(Faxa_swvdf,1),1,6)', ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
end
k=plot(sw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)
title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);

axes(ha(2));
switch n
    case 'CAM6'
        j=plot(repelem(mean(Faxa_rainl,1),1,3)'+ ...
            repelem(mean(Faxa_rainc,1),1,3)'+ ...
            repelem(mean(Faxa_snowl,1),1,3)'+ ...
            repelem(mean(Faxa_snowc,1),1,3)', ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
    case 'CAM4'
        j=plot(repelem(mean(Faxa_rainl,1),1,6)'+ ...
            repelem(mean(Faxa_rainc,1),1,6)'+ ...
            repelem(mean(Faxa_snowl,1),1,6)'+ ...
            repelem(mean(Faxa_snowc,1),1,6)', ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
end
k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)

axes(ha(3));
switch n
    case 'CAM6'
        j=plot(repelem(mean(Sa_tbot,1),1,3), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
    case 'CAM4'
        j=plot(repelem(mean(Sa_tbot,1),1,6), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
end
k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

ylabel('Temperature (K)','FontSize', fontpt)

axes(ha(4));
switch n
    case 'CAM6'
        j=plot(repelem(mean(Sa_shum,1),1,3), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
    case 'CAM4'
        j=plot(repelem(mean(Sa_shum,1),1,6), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
end
k=plot(q_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

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
switch n
    case 'CAM6'
        j=plot(mean(Sa_u,1), ...
            '--', 'color',rgb('darksalmon'), 'linewidth', mean_width); hold on
        k=plot(mean(Sa_v,1), ...
            '-', 'color',rgb('darksalmon'), 'linewidth', mean_width);
        l=plot( (mean(Sa_u,1).^2+ ...
            mean(Sa_v,1).^2).^0.5, ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width);
    case 'CAM4'
        j=plot(repelem(mean(Sa_u,1),1,6)', ...
            '--', 'color',rgb('darksalmon'), 'linewidth', mean_width); hold on
        k=plot(repelem(mean(Sa_v,1),1,6)', ...
            '-', 'color',rgb('darksalmon'), 'linewidth', mean_width);
        l=plot((repelem(mean(Sa_u,1),1,6)'.^2+ ...
            repelem(mean(Sa_v,1),1,6)'.^2).^0.5, ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width);
end
m=plot(wind_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k,l,m], ['Zonal ',n,' Ens Mean'], ['Meridional ',n, ' Ens Mean'], ...
    ['Total ',n,' Ens Mean'],'Tower Met','FontSize',fontpt);

ylabel('Wind (m s^{-1})','FontSize', fontpt)
title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);


axes(ha(2));
switch n
    case 'CAM6'
        j=plot(repelem(mean(Sa_tbot,1),1,3), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
    case 'CAM4'
        j=plot(repelem(mean(Sa_tbot,1),1,6), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
end
k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

ylabel('Temperature (K)','FontSize', fontpt)

axes(ha(3));
j=plot(repelem(mean(Faxa_rainl,1),1,3)'+ ...
    repelem(mean(Faxa_rainc,1),1,3)'+ ...
    repelem(mean(Faxa_snowl,1),1,3)'+ ...
    repelem(mean(Faxa_snowc,1),1,3)', ...
    '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

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
switch n
    case 'CAM6'
        j=plot(repelem(mean(Sa_pbot,1),1,3), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
    case 'CAM4'
        j=plot(repelem(mean(Sa_pbot,1),1,6), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
end
k=plot(ps_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

ylabel('Pressure (Pa)','FontSize', fontpt)
title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);

axes(ha(2));
switch n
    case 'CAM6'
        j=plot(repelem(mean(Faxa_lwdn,1),1,3), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
    case 'CAM4'
        j=plot(repelem(mean(Faxa_lwdn,1),1,6), ...
            '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
end
j=plot(repelem(mean(Faxa_lwdn,1),1,3), ...
    '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
k=plot(lw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);

ylabel('LW radiation (W m^{-1})','FontSize', fontpt)
linkaxes(ha(1:2), 'x');

display('  ')
display(['Finished tower ',n,' diagnostic plots. Press enter to proceed to adjusted ',n,' plots..'])
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
switch n
    case 'CAM6'
        ax1=plot(Faxa_swndr(:,188*24:195*24)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Faxa_swndr_adjust(:,188*24:195*24)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        j=plot(mean(Faxa_swndr(:,188*24:195*24),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        k=plot(mean(Faxa_swndr_adjust(:,188*24:195*24),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        ax1=plot(Faxa_swndr(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Faxa_swndr_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        j=plot(mean(Faxa_swndr(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        k=plot(mean(Faxa_swndr_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
legend( [j,k], [n,' Mean'], ['Adjusted ',n,' Mean'],'FontSize',fontpt);
ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
title([n,' scaling against tower NR1 Met: July ', yearstr] , 'FontSize',fontpt);
switch n
    case 'CAM6'
        xticks([24:24:24*18]);
    case 'CAM4'
end
set(ha(1),'xticklabel',[])

axes(ha(2));
switch n
    case 'CAM6'
        ax2=plot(repelem(Faxa_rain(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Faxa_rainl_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Faxa_rain(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Faxa_rainl_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        ax2=plot(Faxa_rain(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Faxa_rainl_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Faxa_rain(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Faxa_rainl_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Rain (mm  s^{-1})','FontSize', fontpt)
switch n
    case 'CAM6'
        xticks([24:24:24*18]);
    case 'CAM4'
end
set(ha(2),'xticklabel',[])

axes(ha(3));
switch n
    case 'CAM6'
        ax3=plot(repelem(Sa_tbot(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Sa_tbot_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Sa_tbot(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Sa_tbot_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4';
        ax3=plot(Sa_tbot(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Sa_tbot_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Sa_tbot(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Sa_tbot_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Temperature (K)','FontSize', fontpt)
switch n
    case 'CAM6'
        xticks([24:24:24*18]);
    case 'CAM4'
end
set(ha(3),'xticklabel',[])

axes(ha(4));
switch n
    case 'CAM6'
        ax4=plot(repelem(Sa_shum(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Sa_shum_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Sa_shum(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Sa_shum_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        ax4=plot(Sa_shum(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Sa_shum_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Sa_shum(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Sa_shum_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
switch n
    case 'CAM6'
        xticks([24:24:24*18]);
    case 'CAM4'
end
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
switch n
    case 'CAM6'
        plot(Faxa_swndr(:,32*24:45*24)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Faxa_swndr_adjust(:,32*24:45*24)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        j=plot(mean(Faxa_swndr(:,32*24:45*24),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        k=plot(mean(Faxa_swndr_adjust(:,32*24:45*24),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        plot(Faxa_swndr(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Faxa_swndr_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        j=plot(mean(Faxa_swndr(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        k=plot(mean(Faxa_swndr_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
legend( [j,k], [n,' Mean'], ['Adjusted ',n,'  Mean'],'FontSize',fontpt);
ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
title([n,' scaling against tower NR1 Met: February ', yearstr] , 'FontSize',fontpt);
switch n
    case 'CAM6'
        xticks([24:24:24*13]);
    case 'CAM4'
        xticks([6:6:6*13]);
end

axes(ha(2));
switch n
    case 'CAM6'
        plot(repelem(Faxa_snow(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Faxa_snowl_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Faxa_snow(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Faxa_snowl_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        plot(Faxa_snow(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Faxa_snowl_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Faxa_snow(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Faxa_snowl_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Snow (mm  s^{-1})','FontSize', fontpt)
switch n
    case 'CAM6'
        xticks([24:24:24*13]);
    case 'CAM4'
        xticks([6:6:6*13]);
end

axes(ha(3));
switch n
    case 'CAM6'
        plot(repelem(Sa_tbot(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Sa_tbot_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Sa_tbot(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Sa_tbot_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        plot(Sa_tbot(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Sa_tbot_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Sa_tbot(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Sa_tbot_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Temperature (K)','FontSize', fontpt)
switch n
    case 'CAM6'
        xticks([24:24:24*13]);
    case 'CAM4'
        xticks([6:6:6*13]);
end

axes(ha(4));
switch n
    case 'CAM6'
        plot(repelem(Sa_shum(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Sa_shum_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Sa_shum(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Sa_shum_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        plot(Sa_shum(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Sa_shum_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Sa_shum(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Sa_shum_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
switch n
    case 'CAM6'
        xticks([24:24:24*13]);
    case 'CAM4';
        xticks([6:6:6*13]);
end
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
switch n
    case 'CAM6'
        k=plot(sw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
    case 'CAM4'
        k=plot(sw_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
end
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);
title([n,' scaling against tower NR1 Met ', yearstr] , 'FontSize',fontpt);
ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)


axes(ha(2));

ax1=plot(Faxa_swndr', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
plot(Faxa_swndr_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
j=plot(mean(Faxa_swndr,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
k=plot(mean(Faxa_swndr_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
grid on
legend( [j,k], [n,' Mean'], [n,' Adjusted Mean'],'FontSize',fontpt);
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
switch n
    case 'CAM6'
        k=plot(ppt_3hr, '-', 'color',rgb('black'),'linewidth', mean_width);
    case 'CAM4'
        k=plot(ppt_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
end
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);
title([n,' scaling against tower NR1 Met ', yearstr] , 'FontSize',fontpt);
ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)
ylim([0 2*10^-3])

axes(ha(2));
j=plot(mean(Sa_tbot,1), ...
    '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
switch n
    case 'CAM6'
        k=plot(ta_3hr, '-', 'color',rgb('black'),'linewidth', mean_width);
    case 'CAM4'
        k=plot(ta_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
end
grid on
legend([j,k], [n,' Ensemble Mean'], 'Tower Met','FontSize',fontpt);
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
switch n
    case 'CAM6'
        m=plot(wind_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
    case 'CAM4'
        m=plot(wind_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
end
grid on
legend([j,k,l,m], ['Zonal ',n,' Ens Mean'], ['Meridional ',n,' Ens Mean'], ...
    ['Total, ',n,' Ens Mean'],'Tower Met','FontSize',fontpt);

ylabel('Wind (m s^{-1})','FontSize', fontpt)
title(['Niwot Ridge Meteorology ', yearstr] , 'FontSize',fontpt);
axes(ha(2));

ax1=plot(Sa_u', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
plot(Sa_u_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);

j=plot(mean(Sa_u,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
k=plot(mean(Sa_u_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
grid on
legend( [j,k], [n,' Mean'], ['Adjusted ',n,' Mean'],'FontSize',fontpt);
ylabel('Zonal wind (m s^{-1})','FontSize', fontpt)
axes(ha(3));

ax2=plot(Sa_v', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
plot(Sa_v_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
plot(mean(Sa_v,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
plot(mean(Sa_v_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
grid on
ylabel('Meridional wind (m s^{-1})','FontSize', fontpt)
axes(ha(4));
switch n
    case 'CAM6'
        ax3=plot(repelem(Sa_pbot,1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Sa_pbot_adjust,1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Sa_pbot,1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Sa_pbot_adjust,1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        ax3=plot(Sa_pbot', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Sa_pbot_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Sa_pbot,1), '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Sa_pbot_adjust,1), '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Pressure (Pa)','FontSize', fontpt)

axes(ha(5));
switch n
    case 'CAM6'
        ax4=plot(repelem(Faxa_lwdn,1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(repelem(Faxa_lwdn_adjust,1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(repelem(mean(Faxa_lwdn,1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(repelem(mean(Faxa_lwdn_adjust,1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
    case 'CAM4'
        ax4=plot(Faxa_lwdn', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
        plot(Faxa_lwdn_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
        plot(mean(Faxa_lwdn,1), '-', 'color',rgb('darkred'), 'linewidth', mean_width);
        plot(mean(Faxa_lwdn_adjust,1), '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
end
grid on
ylabel('Longwave Radiation (W m^{-2})','FontSize', fontpt)
linkaxes(ha(1:5), 'x');

display('   ')
display(['Adjusted ',n,' plots and diagnostics completed.  Press enter to write ',yearstr,' ',n,' adjusted files to netcdf'])
display('   ')
pause
close all


end
