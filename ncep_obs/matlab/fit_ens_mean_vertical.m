% obsfit.m plot of the error of analysis and guess
clf
p_v=load('Tges_ver_ave.dat');
a_v=load('Tanl_ver_ave.dat');
yp_v=p_v(:,1);
ya_v=a_v(:,1);
xp_v=p_v(:,2);
xa_v=a_v(:,2);
%
subplot('position', [0.1,0.6,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T rmse, ensemble mean, NH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v=p_v(:,4);
xa_v=a_v(:,4);
%
subplot('position', [0.6,0.6,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T rmse, ensemble mean, SH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v=p_v(:,6);
xa_v=a_v(:,6);
subplot('position', [0.1,0.1,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T rmse, ensemble mean, TR', 'fontsize', 10)
%
xp_v=p_v(:,8);
xa_v=a_v(:,8);
subplot('position', [0.6,0.1,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T rmse, ensemble mean, NA', 'fontsize', 10)
%
print -dpsc t_vertical.ps
pause
%
%
% - ----------------- W --------------
p_v=load('Wges_ver_ave.dat');
a_v=load('Wanl_ver_ave.dat');
yp_v=p_v(:,1);
ya_v=a_v(:,1);
xp_v=p_v(:,2);
xa_v=a_v(:,2);
%
subplot('position', [0.1,0.6,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('Vector wind rms error, ensemble mean, NH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v=p_v(:,4);
xa_v=a_v(:,4);
%
subplot('position', [0.6,0.6,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('Vector wind rmse, ensemble mean, SH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v=p_v(:,6);
xa_v=a_v(:,6);
subplot('position', [0.1,0.1,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('Vector wind rmse, ensemble mean, TR', 'fontsize', 10)
%
xp_v=p_v(:,8);
xa_v=a_v(:,8);
subplot('position', [0.6,0.1,0.35,0.35]), plot(xp_v,yp_v,'c+-',xa_v,ya_v,'ro:')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('Vector wind rmse, ensemble mean, NA', 'fontsize', 10)
%
print -dpsc w_vertical.ps
