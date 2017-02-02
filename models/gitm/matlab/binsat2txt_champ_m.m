%% CHAMP-flown-through-GITM data reader (from GITM binary)

% DART $Id$
% CREDIT: Alexey Morozov

%run champ2txt.m first

% clear
% clc
% close all

format compact

% cd ~/Downloads/sats/data
% cd /Users/admin/gitm/GITM2/run/UA/data/
% cd /Applications/HEAD/DART9/models/gitm/work/advance_temp_e4/UA/data
% cd /Applications/HEAD/DART10/models/gitm/work/advance_temp_e4/UA/data



%  ns=1; %NUMBER OF SATELLITES
%  sns=num2str(ns); %name of the output text file is ['text' sns '.txt']
idt=fopen('text.txt','w');

[~,wlon]=min(abs(LongitudeT-(96.2626+180) )) %lon closest to AA (WI in 9x9)
[~,wlat]=min(abs(LatitudeT-42.2526)) %lat closest to AA (WI in 9x9)

% wlon=5;
% wlat=14;
walt=35; 

LongitudeT(wlon)
LatitudeT(wlat)
AltitudeT(walt)


LonCi=interp1(tc,LongitudeChamp2,tg);
LatCi=interp1(tc,LatitudeChamp2,tg);
AltCi=interp1(tc,AltitudeChamp2,tg);


Alt38=nan(size(tg));
AltitudeG2=nan(size(tg));
RhoG38=nan(size(tg));
RhoG2=nan(size(tg));
LongitudeSat=nan(size(tg));
LatitudeSat=nan(size(tg));
AltitudeSat=nan(1,nAlts);
RhoGITM=nan(1,nAlts);
Rw=nan(size(tg));

yy=nan(size(tg));
mo=yy;
dd=yy;
hh=yy;
mm=yy;
ss=yy;

for m=1:length(tg) %iterate over number of minutes
    [yy(m) mo(m) dd(m) hh(m) mm(m) ss(m)]=datevec(datenum(2000+y2,0,0)+d1+tg(m)/1440);
    
    
    fname=['cham_001_t' ...
        num2str(yy(m)-2000,'%02.0f') ...
        num2str(mo(m),'%02.0f') ...
        num2str(dd(m),'%02.0f') ...
        '_' ...
        num2str(hh(m),'%02.0f') ...
        num2str(mm(m),'%02.0f') ...
        num2str(ss(m),'%02.0f') ...
        '.sat'];
    
    disp(fname)
    
    idf=fopen(fname,'r');
    
    
    
    for iAlt=1:54 %iterate over altitudes
        nvars=fread(idf,1,'int32')/8;
        temp=fread(idf,nvars,'float64');
        nvars=fread(idf,1,'int32')/8;
        
        AltitudeSat(iAlt)     = temp(3); %km
        RhoGITM(iAlt)          = temp(4);
    end
    fclose(idf);
    
    % disp(        AltitudeSat(38))
    % disp(' ')
    Alt38(m) = AltitudeSat(38);
    RhoG38(m) = RhoGITM(38);
    
    %%simple way (no interpolation, just closest altitude/time)
    %        [~, sec_ind] = min(abs(sec(:) - m(1)*60));
    %        [~, a_ind  ] = min(abs(Altitude(3:52) - AltitudeChamp(sec_ind) ));
    %        a_ind=a_ind+2; %add 2 because search array is an exerpt
    %        AltitudeG1(m)=Altitude(a_ind);
    %        RhoG1(m)=RhoGITM(a_ind);
    
    %%harder - interpolate altitude
    AltitudeG2(m) = AltCi(m);
    RhoG2(m)=interp1(AltitudeSat,RhoGITM,AltCi(m));
    
    %%hardest - interpolate altitude AND then time %need to search over m and
    %%save Alt(m,ialt). meh, don't feel like it
    %         [~, sec_ind] = min(abs(sec(:) - m(1)*60)); %m(1): (1) is to remind me that m is a scalar and sec(:) is the array I'm searching over
    %         a_ind = find(Altitude(3:52) < AltitudeChamp(sec_ind) );
    %         a_ind = a_ind(end) + 2; % last one smaller than the champ and add 2 because search array is an exerpt
    %         a = ( Altitude(a_ind+1) - Altitude(a_ind) )/( Altitude(a_ind+1)-Altitude(a_ind) ); %I know this is always 1, just a (in?)sanity check
    %         b = ( Altitude(a_ind)*Altitude(a_ind+1) - Altitude(a_ind+1)*Altitude(a_ind) )/( Altitude(a_ind+1)-Altitude(a_ind) ); %zero
    %         AltitudeG2(m) = a*AltitudeChamp(sec_ind) + b ;
    %         a = ( RhoGITM(a_ind+1) - RhoGITM(a_ind) )/( Altitude(a_ind+1)-Altitude(a_ind) ) ;
    %         b = ( RhoGITM(a_ind)*Altitude(a_ind+1) - RhoGITM(a_ind+1)*Altitude(a_ind) )/( Altitude(a_ind+1)-Altitude(a_ind) ) ;
    %         RhoG3(m) = a*AltitudeChamp(sec_ind) + b ;
    
    %%hardester - interpolate altitude and time %implausible in this loop as
    %%Alt is not stored over time.
    %         si = find(m(:) < sec(ii)/60 );
    %         ai = find(Altitude(3:52) < AltitudeChamp(jj) );
    %         ai = ai(end) + 2;
    %         abc=[Altitude(ai:ai+2)' (m(si:si+2))' ones(3,1)]\RhoGITM(ai:ai+2)';
    %         RhoG4(m) = abc*[AltitudeChamp(jj); sec(ii)/60 ;  1];
    
    LongitudeSat(m)    = temp(1)*180/pi; %for these numbers see  3DALL_... .header file
    LatitudeSat(m)     = temp(2)*180/pi;
    %         Altitude(m)     = temp(3);
    %        Temperature(m)  = temp(16);
    %        NDSo3(m)        = temp(5);
    %        IDSe(m)         = temp(34);
    %         RhoGITM(m)          = temp(4);
    
    
    
    
    Rw(m)=RhoT(wlon,wlat,walt,m);
    
    
end

%% .txt for ./t2o iChamp loc

sp=ones(length(tg),1);


text_M=[num2str(kind*sp) 32*sp ...
    num2str(LatitudeSat','%9.5f') 32*sp ...
    num2str(LongitudeSat','%10.5f') 32*sp ...
    num2str(AltitudeG2','%8.1f') 32*sp ...
    num2str(yy') 32*sp ...
    num2str(mo','%02.0f') 32*sp ...
    num2str(dd','%02.0f') 32*sp ...
    num2str(hh','%02.0f') 32*sp ...
    num2str(mm','%02.0f') 32*sp ...
    num2str(ss','%02.0f') 32*sp ...
    num2str(RhoG2','%12.6e') 32*sp ...
    num2str(v*sp) 32*sp ...
    10*sp];

fwrite(idt,text_M');

%% .txt for ./t2o 1 loc

% sp=ones(length(tg),1);
% 
% 
% text_M=[num2str(kind*sp) 32*sp ...
%     num2str(LatitudeT(wlat)*sp,'%9.5f') 32*sp ...
%     num2str(LongitudeT(wlon)*sp,'%10.5f') 32*sp ...
%     num2str(AltitudeT(walt)*sp,'%8.1f') 32*sp ...
%     num2str(yy') 32*sp ...
%     num2str(mo','%02.0f') 32*sp ...
%     num2str(dd','%02.0f') 32*sp ...
%     num2str(hh','%02.0f') 32*sp ...
%     num2str(mm','%02.0f') 32*sp ...
%     num2str(ss','%02.0f') 32*sp ...
%     num2str(Rw','%12.6e') 32*sp ...
%     num2str(v*sp) 32*sp ...
%     10*sp];
% 
% fwrite(idt,text_M');


%% .txt for ./t2o 1 fancy loc % if you want to write obs_seq.out at some fancy location (like subsolar point)

%Ann Arbor
% LonAA=277.5 *ones(size(tg)); %82.5d W, need to repeat it tg times, because it is defined over tg and doesn't change over tg
% LatAA=42.5  *ones(size(tg)); %42.5d N
% AltAA=393983.5 *ones(size(tg)); %about 400km
% 
% RhoAA=interpolateN( RhoT(:,:,:,:), ... %doesn't work too well if you're interpolating to the gridpoint (it takes the cube closest to the gridpoint, if you want perfection at gridpoints, use interpolateN_b)
%     {LongitudeT;LatitudeT;AltitudeT;tg}, ...
%     {LonAA;LatAA;AltAA;tg}); 

%Subsolar Point (or just equator)
LonAA=mod( (1440-mod(tg,1440))/1440*360+180 , 360); %subsolar point (point on Earth closest to the Sun) moves Westward at 360deg/day=15deg/hr. It is where the local noon is (ie where it is 12:00 military time, not 00:00).
LatAA=  (23+26/60)*sin(2*pi/365*((d1+tg/1440)-78) ); %subsolar point moves between the Tropic of Cancer (Northern Solstice, let's say Jun 22) and the Tropic of Capricorn (Southern Solstice, let's say Dec 21). 
%+ so it is a sinusoid with magnitude of +-23.43deg (that's where the Tropics are) and frequency of 2pi rev/365 days. The phase shift is approximately datenum(2002,3,20)-datenum(2002,1,1) 
AltAA=393983.5 *ones(size(tg)); %about 400km

RhoAA=interpolateN( RhoT(:,:,:,:), ... %doesn't work too well if you're interpolating to the gridpoint (it takes the cube closest to the gridpoint, if you want perfection at gridpoints, use interpolateN_b)
    {LongitudeT;LatitudeT;AltitudeT;tg}, ...
    {LonAA;LatAA;AltAA;tg}); 


% plot(tg,RhoAA,'r')
% hold on
% plot(tg,squeeze(RhoT(wlon,wlat,walt,:)), 'b')
% hold off
% ylim([0 10e-12])


% sp=ones(length(tg),1);
% 
% 
% text_M=[num2str(kind*sp) 32*sp ... %for txt2obs
%     num2str(LatAA','%9.5f') 32*sp ...
%     num2str(LonAA','%10.5f') 32*sp ...
%     num2str(AltAA','%8.1f') 32*sp ...
%     num2str(yy') 32*sp ...
%     num2str(mo','%02.0f') 32*sp ...
%     num2str(dd','%02.0f') 32*sp ...
%     num2str(hh','%02.0f') 32*sp ...
%     num2str(mm','%02.0f') 32*sp ...
%     num2str(ss','%02.0f') 32*sp ...
%     num2str(RhoAA','%12.6e') 32*sp ...
%     num2str(v*sp) 32*sp ...
%     10*sp];
%     
%     fwrite(idt,text_M');

% text_M=[num2str(yy') 32*sp ... %for gitm sat.dat
%     num2str(mo','%02.0f') 32*sp ...
%     num2str(dd','%02.0f') 32*sp ...
%     num2str(hh','%02.0f') 32*sp ...
%     num2str(mm','%02.0f') 32*sp ...
%     num2str(ss','%02.0f') 32*sp ...
%     48*sp 32*sp ...
%     num2str(LonAA','%10.5f') 32*sp ...
%     num2str(LatAA','%9.5f') 32*sp ...
%     num2str(AltAA'/1000,'%7.3f') 32*sp ...    
%     num2str(RhoAA','%12.6e') 32*sp ...
%     10*sp];
% fwrite(idt,['File made WHEN  : ' datestr(now) ' ' 10]);
% fwrite(idt,['File made HOW   : matlab subsolar point calculation for days' num2str(d1) '+' num2str(nd) ' ' 10]);
% fwrite(idt,['File made WHERE : ' pwd  10]);
% fwrite(idt,['Matlab code : binsat2txt_champ_m.m '  10]);
% fwrite(idt,['#START' 10]);
% 
% fwrite(idt,text_M');

%%
fclose(idt);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
