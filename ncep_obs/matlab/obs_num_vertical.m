% obsfit.m plot of the error of analysis and guess
clf
p_v=load('Tges_ver_ave.dat');
yp_v=p_v(:,1);
xp_v_num=p_v(:,3);
%
subplot('position', [0.1,0.6,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T num, NH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v_num=p_v(:,5);
%
subplot('position', [0.6,0.6,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T num, SH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v_num=p_v(:,7);
subplot('position', [0.1,0.1,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T num, TR', 'fontsize', 10)
%
xp_v_num=p_v(:,9);
subplot('position', [0.6,0.1,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%axis([0 10 100 1000])
grid
set(gca,'YDir', 'reverse')
%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('T num, NA', 'fontsize', 10)
%
print -dpsc t_num_vert.ps
pause
%
%
% - ----------------- W --------------
p_v=load('Wges_ver_ave.dat');
yp_v=p_v(:,1);
xp_v_num=p_v(:,3);
%
subplot('position', [0.1,0.6,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel('wind num,  NH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v_num=p_v(:,5);
%
subplot('position', [0.6,0.6,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel(' wind num, SH', 'fontsize', 10)
legend('guess', 'analysis')
h=legend;
legend(h,'boxoff')
%
xp_v_num=p_v(:,7);
subplot('position', [0.1,0.1,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel(' wind num, TR', 'fontsize', 10)
%
xp_v_num=p_v(:,9);
subplot('position', [0.6,0.1,0.35,0.35]), plot(xp_v_num ,yp_v,'c+-')
%axis([2 15 100 1000])
grid
set(gca,'YDir', 'reverse')
%%ylabel('Pressure(hPa)', 'fontsize', 10)
xlabel(' wind num, NA', 'fontsize', 10)
%
print -dpsc w_num_vert.ps
