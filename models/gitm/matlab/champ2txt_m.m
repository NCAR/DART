%% CHAMP ascii files reader, text.txt writer, and GITM champ_*.dat files writer

% DART $Id$
% CREDIT: Alexey Morozov

tc=[];
LongitudeChamp2s=[];
LatitudeChamp2s=[];
AltitudeChamp2s=[];
RhoChamp2s=[];
v2s=[];

% fn2=['champ_with_rho_Dec1to10.dat']
% idt=fopen(fn2,'w'); %dd(i)
% 
% 
% fwrite(idt,[' '  10]);
% 
% fwrite(idt,['File made on : ' datestr(now) ' ' 10]);
% fwrite(idt,['File made FROM : Density_3deg_02_335+ ' 10]);
% fwrite(idt,['on sisko.colorado.edu/sutton/data/ver2.2/champ/density/2002/ascii/ ' 10]);
% fwrite(idt,['Matlab code : champ2txt.m '  10]);
% 
% fwrite(idt,[' '  10]);
% fwrite(idt,['#START '  10]);

%% .txt for ./t2o
idt=fopen('text.txt','w');

for d=1:nd
    
    fn=['Density_3deg_' num2str(y2,'%02.0f') '_' num2str(d1-1+d,'%03.0f') '.ascii']
    % fn='Density_3deg_02_335.ascii';
    
    [s,r] = unix(['wc -l ' fn]);
    
    r=str2num(r(1:find(r>57,1)-1)); %take only the number part of reply
    rs(d)=r;
    unix(['tail -n ' num2str(r-2) ' ' fn ' > D_trunk.dat']);
    s=load('D_trunk.dat');
%     s(1,:)
    
    % Version 2.2;
    % 	Ascii files are arranged in columnn format with a 1-line header of Descriptions
    % 	including units used) separated by semi-colons (;). The length of the files is
    % 	approximately 1840 lines.
    %
    % 	Column:			Format:		Description:		Units:
    % 	--------------------------------------------------------------------------------
    % 	[1-2]			2I      1   Two-digit Year		years
    % 	[4-6]			3I      2   Day of the Year		days
    % 	[8-15]			8.3F	3	Second of the Day	GPS time,sec
    % 	[17-19]			3I      4   Latitude Center of 	deg Current Averaging Bin
    % 	[21-29]			9.5F	5	Geodetic Latitude	deg
    % 	[31-40]			10.5F	6	Longitude		deg
    % 	[42-48]			7.3F	7	Satellite Height	km
    % 	[50-56]			7.4F	8	Satellite Local Time	hours
    % 	[58-66]			9.5F	9	Quasi-Dipole Latitude	deg
    % 	[68-77]			10.5F	10	Magnetic Longitude	deg
    % 	[79-85]			7.4F	11	Magnetic Local Time	hours
    % 	[87-98]			12.6E	12	Neutral Density		kg/m^3
    % 	[100-111]		12.6E	13	Density	@ 400km		kg/m^3
    % 	[113-124]		12.6E	14	Density @ 410km		kg/m^3
    % 	[126-137]		12.6E	15	NRL-MSIS Density at	kg/m^3 Satellite Height
    % 	[139-150]		12.6E	16	Density Uncertainty	kg/m^3
    % 	[152-153]		2I		17  Number of Data Points in Current Averaging Bin
    % 	[155-156]		2I		18  Number of Points that Require Interpolation
    % 	[158-162]       5.3F	19	Coefficient of Drag Averaged over Bin
    
    [yy mo dd hh mm ss]=datevec(datenum(2000+s(:,1),0,0)+s(:,2)+s(:,3)/86400);
    disp(datevec(datenum(2000+s(1,1),0,0)+s(1,2)+s(1,3)/86400))
    ss=round(ss);% round because fractional seconds are not supported in DART
    
    sec2=s(:,3);
    % for i=1:length(sec2)
    %     %     if i<=r/2 % hack for 2 day lump file
    %     %         dd(i)=20;
    %     %     else
    %     %         dd(i)=21;
    %     %     end
    %     dd(i)=21; %%%%%%%%HHHHHHHHHAAAARD CODED
    %
    %     hh(i)=floor(sec2(i)/3600);
    %     mm(i)=floor(mod(sec2(i),3600)/60);
    %     ss(i)=mod(sec2(i),60);
    % end
    
    % tp = find(hh==12); %time period 917:992
    tp = 1:r-2;        %all
    % tp = 1:380; %0-4:45am
    % tp = 1:78; %0-1am
    
    yy=yy(tp);
    mo=mo(tp);
    dd=dd(tp);
    hh=hh(tp);
    mm=mm(tp);
    ss=ss(tp);
    sec2=sec2(tp);
    
    LongitudeChamp2    = s(tp,6); %lon -180 180
    LatitudeChamp2     = s(tp,5);
    
    AltitudeChamp2     = s(tp,7)*1000;
    RhoChamp2          = s(tp,12);
    RhoMsis2      = s(tp,15);
    v2            = s(tp,16);
    
    
    test=LongitudeChamp2<0;
    LongitudeChamp2(test)=LongitudeChamp2(test)+360;
    
        tc=[tc sec2'+(d-1)*86400];
        LongitudeChamp2s=[LongitudeChamp2s; LongitudeChamp2];
        LatitudeChamp2s=[LatitudeChamp2s; LatitudeChamp2];
        AltitudeChamp2s=[AltitudeChamp2s; AltitudeChamp2];
        RhoChamp2s=[RhoChamp2s; RhoChamp2];
        v2s=[v2s; v2];
    

    
    %% .txt for ./t2o

    sp=ones(r-2,1);    

    
    text_M=[num2str(kind*sp) 32*sp ...
        num2str(LatitudeChamp2,'%9.5f') 32*sp ...
        num2str(LongitudeChamp2,'%10.5f') 32*sp ...
        num2str(AltitudeChamp2,'%8.1f') 32*sp ... %meters !!
        num2str(yy) 32*sp ... 
        num2str(mo,'%02.0f') 32*sp ...
        num2str(dd,'%02.0f') 32*sp ...
        num2str(hh,'%02.0f') 32*sp ... 
        num2str(mm,'%02.0f') 32*sp ...
        num2str(ss,'%02.0f') 32*sp ...
        num2str(RhoChamp2,'%12.6e') 32*sp ...
        num2str(v2) 32*sp  ...
        10*sp];
    
    fwrite(idt,text_M');
    
end

fclose(idt);

LongitudeChamp2=LongitudeChamp2s;
LatitudeChamp2 =LatitudeChamp2s;
AltitudeChamp2 =AltitudeChamp2s;
RhoChamp2 =RhoChamp2s;
v2 =v2s;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
