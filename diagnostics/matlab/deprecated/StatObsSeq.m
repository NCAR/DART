%% StatObsSeq
%
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

QTY_U = 1;
QTY_V = 2;
QTY_T = 4;

a = ReadObsSeq('obs_seq.final');

num_true_times = 21;

iu  = zeros(1,2*num_true_times);
iv  = zeros(1,2*num_true_times);
it  = zeros(1,2*num_true_times);
bu  = zeros(1,2*num_true_times);
bv  = zeros(1,2*num_true_times);
bt  = zeros(1,2*num_true_times);
ru  = zeros(1,2*num_true_times);
rv  = zeros(1,2*num_true_times);
rt  = zeros(1,2*num_true_times);
x   = zeros(1,2*num_true_times);
omu = zeros(1,2*num_true_times);
omv = zeros(1,2*num_true_times);
omt = zeros(1,2*num_true_times);
tmu = zeros(1,2*num_true_times);
tmv = zeros(1,2*num_true_times);
tmt = zeros(1,2*num_true_times);
su  = zeros(1,2*num_true_times);
sv  = zeros(1,2*num_true_times);
st  = zeros(1,2*num_true_times);

days = -1;
secs = -1;
itime = 0;

for i = 1:a.num_obs,

   if (a.days(i) ~= days | a.secs(i) ~= secs)
     itime = itime + 1;
     days = a.days(i);
     secs = a.secs(i);
     x(2*itime-1) = days*86400 + secs;
     x(2*itime)   = days*86400 + secs;
   end

   if ( a.kind(i) == QTY_U )

     if (a.obs(1,i) == -888888.0)

        bu(2*itime-1) = bu(2*itime-1) + 1;
        bu(2*itime)   = bu(2*itime)   + 1;

     elseif (a.qc(1,i) ~= 0.0)

        ru(2*itime-1) = ru(2*itime-1) + 1;
        ru(2*itime)   = ru(2*itime)   + 1;

     else

         iu(2*itime-1) =  iu(2*itime-1) + 1;
         iu(2*itime)   =  iu(2*itime)   + 1;
        omu(2*itime-1) = omu(2*itime-1) + (a.obs(1,i) - a.obs(3,i)).^2;
        omu(2*itime)   = omu(2*itime)   + (a.obs(1,i) - a.obs(4,i)).^2;
        tmu(2*itime-1) = tmu(2*itime-1) + (a.obs(2,i) - a.obs(3,i)).^2;
        tmu(2*itime)   = tmu(2*itime)   + (a.obs(2,i) - a.obs(4,i)).^2;
         su(2*itime-1) =  su(2*itime-1) +  a.obs(5,i).^2;
         su(2*itime)   =  su(2*itime)   +  a.obs(6,i).^2;

     end

   elseif ( a.kind(i) == QTY_V )

     if (a.obs(1,i) == -888888.0)

        bv(2*itime-1) = bv(2*itime-1) + 1;
        bv(2*itime)   = bv(2*itime)   + 1;

     elseif (a.qc(1,i) ~= 0.0)

        rv(2*itime-1) = rv(2*itime-1) + 1;
        rv(2*itime)   = rv(2*itime)   + 1;

     else

         iv(2*itime-1) =  iv(2*itime-1) + 1;
         iv(2*itime  ) =  iv(2*itime  ) + 1;
        omv(2*itime-1) = omv(2*itime-1) + (a.obs(1,i) - a.obs(3,i)).^2;
        omv(2*itime  ) = omv(2*itime  ) + (a.obs(1,i) - a.obs(4,i)).^2;
        tmv(2*itime-1) = tmv(2*itime-1) + (a.obs(2,i) - a.obs(3,i)).^2;
        tmv(2*itime  ) = tmv(2*itime  ) + (a.obs(2,i) - a.obs(4,i)).^2;
         sv(2*itime-1) =  sv(2*itime-1) +  a.obs(5,i).^2;
         sv(2*itime  ) =  sv(2*itime  ) +  a.obs(6,i).^2;

     end

   elseif ( a.kind(i) == QTY_T )

     if (a.obs(1,i) == -888888.0)

        bt(2*itime-1) = bt(2*itime-1) + 1;
        bt(2*itime)   = bt(2*itime  ) + 1;

     elseif (a.qc(1,i) ~= 0.0)

        rt(2*itime-1) = rt(2*itime-1) + 1;
        rt(2*itime)   = rt(2*itime  ) + 1;

     else

         it(2*itime-1) =  it(2*itime-1) + 1;
         it(2*itime  ) =  it(2*itime  ) + 1;
        omt(2*itime-1) = omt(2*itime-1) + (a.obs(1,i) - a.obs(3,i)).^2;
        omt(2*itime  ) = omt(2*itime  ) + (a.obs(1,i) - a.obs(4,i)).^2;
        tmt(2*itime-1) = tmt(2*itime-1) + (a.obs(2,i) - a.obs(3,i)).^2;
        tmt(2*itime  ) = tmt(2*itime  ) + (a.obs(2,i) - a.obs(4,i)).^2;
         st(2*itime-1) =  st(2*itime-1) +  a.obs(5,i).^2;
         st(2*itime  ) =  st(2*itime  ) +  a.obs(6,i).^2;

     end

   end

end

num_true_times = itime;

omu = sqrt(omu./iu);
omv = sqrt(omv./iv);
omt = sqrt(omt./it);
tmu = sqrt(tmu./iu);
tmv = sqrt(tmv./iv);
tmt = sqrt(tmt./it);
 su = sqrt( su./iu);
 sv = sqrt( sv./iv);
 st = sqrt( st./it);

x = x - x(1);
time_unit = 'seconds';
if (max(x) > 60.0)
     x = x/60;
     time_unit = 'minutes';
end
if (max(x) > 60.0)
     x = x/60;
     time_unit = 'hours';
end
if (max(x) > 24.0)
     x = x/24;
     time_unit = 'days';
end

subplot(3,1,1);
plot(x,tmu,x,omu,x,su)
title('O-U RMS (m/s)')
xlim([(x(1)-0.5) (max(x)+0.5)])

subplot(3,1,2);
plot(x,tmv,x,omv,x,sv)
title('O-V RMS (m/s)')
xlim([(x(1)-0.5) (max(x)+0.5)])

subplot(3,1,3);
plot(x,tmt,x,omt,x,st)
title('O-T RMS (K)')
xlim([(x(1)-0.5) (max(x)+0.5)])

xlabel(time_unit,'Fontsize',12)

%legend('CV3','CV5 ENS')
legend('RMS diff truth','RMS diff obs','Spread')

figure

subplot(3,1,1);
plot(x,iu,x,iv,x,it)
title('Assimilated Obs')
xlim([(x(1)-0.5) (max(x)+0.5)])

subplot(3,1,2);
plot(x,ru,x,rv,x,rt)
title('Rejected Obs')
xlim([(x(1)-0.5) (max(x)+0.5)])

subplot(3,1,3);
plot(x,bu,x,bv,x,bt)
title('Bad Obs')
xlim([(x(1)-0.5) (max(x)+0.5)])

xlabel(time_unit,'Fontsize',12)

legend('U','V','T')


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

