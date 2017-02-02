function [tts,rs,gs]=loc(tL, LonL, LatL, AltL, ...
                        tV, LonV, LatV, AltV, RhoV, ...
                        cfh, cfv)
%localizing function 
%
% location to which you want to Localize to is defined by LonL (L for LOCATION), 
% LatL, AltL, all defined at common time tL
%
% location where you have Values at is defined by tV, LonV, LatV, AltV, where the value itself is RhoV
%
% Lons and Lats are assumed to be in degrees, Alts in meters
%
% outputs times (tts), values (rs) and localization values (from 0 to 1,
% gs)

% DART $Id$
% CREDIT: Alexey Morozov

dn=abs(LonL-LonV);
dn(dn>180)=360-dn(dn>180); 
dt=LatL-LatV;

phi_s = LatV*pi/180;
lam_s = LonV*pi/180;
phi_f = LatL*pi/180;
lam_f = LonL*pi/180;

ds=atan2(sqrt( ( cos(phi_f).*sin(abs(lam_s-lam_f)) ).^2 + ( cos(phi_s).*sin(phi_f)-sin(phi_s).*cos(phi_f).*cos(abs(lam_s-lam_f)) ).^2) ,( sin(phi_s).*sin(phi_f)+cos(phi_s).*cos(phi_f).*cos(abs(lam_s-lam_f)) ));
% ds=2*asin(sqrt( sin(abs(phi_s-phi_f)/2).^2 + cos(phi_s).*cos(phi_f).*(sin(abs(lam_s-lam_f)/2).^2) ) );
% keyboard

da=( ( AltL-AltV ) /cfv); %/cfv as per Nancy's Oct/5/12 email, 
% da=( ( AltL-AltV ) /cfv)*180/pi; %/cfv as per Nancy's Oct/5/12 email, *180/pi to convert rad to deg (since dn and dt are in deg).

gc=GC(cfh*pi/180,sqrt(ds.^2+da.^2)); %compute 3D straight line distance in degrees and apply Gaspari Cohn localization with cfh cutoff
gc=gc(:); %gc must be a column in order for sort(*,1,*) to work

[gci,j] = sort(gc,1,'ascend'); %sort the result
tci=tV(j); %assumes tL=tV
rci=RhoV(j);

jj=gci>.1; %throw away everything far away
gs=gci(jj); 
tts=tci(jj);
rs=rci(jj); 

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
