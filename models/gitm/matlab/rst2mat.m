% RST2MAT - A program for reading and plotting GITM restart files, which
%   are assumed to be binary and created by fortran90 program running on an
%   Intel Mac. It also assumes that the files are in the current directory
%   and are named "bXXXX.rst", where XXXX goes from 1 to number of blocks.
%
%   It requires file called "imagescn.m" to be in the current directory or in 
%   your Matlab path. IMAGESCN is used to plot the 2D-true-area-scale
%   plots. See "help imagescn".
%
%   This file clears some workspace variables, so see line 30 before running the code.
%   Also, line 35 cd's into your restart folder of interes (change as needed).
%
%   There are 4 sections:
%       1. READING the files
%       2. PUTTING-TOGETHER the data from multiple files
%       3. PLOTTING the data
%       4. (optional) MAT2RST - if you want to write your own restart files for diagnostics or something - read comments in that section
%
%   The variables in the READING section (Lon, Lat, Alt, NDS, Temp, V, etc)
%   have the last dimension of nBlocksMax, whereas the variables in section
%   PUTTING-TOGETHER (LonT, LatT, AltT, NDST, TempT, VT, etc.) don't have 
%   this dimension as the name of section suggests and have a "T" on the
%   end of their name, which stands for Total.
%   The last plot in the PLOTTING section is composed of salvaged pieces of
%   SPHERE3D code written by by J De Freitas:
%   http://www.mathworks.com/matlabcentral/fileexchange/8585-3d-plot-on-a-spherical-surface
%
%   Alexey Morozov 09/01/2011 morozova aatt umich ddott edu

% DART $Id$
% CREDIT: Alexey Morozov

clear %harsh <-
%milder version:
% clear Lon LonT Lat LatT Alt AltT Temp Tempt TempT eTemp eTempt eTempT ITemp ITempt ITempT NDS NDSt NDST IDS IDSt IDST V Vt VT IV IVt IVT VV VVt VVT F107 F107t Rho Rhot RhoT 
clc
close all
format compact
% format long

%% 1.READING section

h0='/Applications/HEAD/DART8/'; %where is the DART folder?
wd='advance_temp_e1/UA/data';     %work directory (ie where the truth run is, if you ran 20 ensemble members, this folder name will end in 21, 3 => 4, n => n+1)
cd([h0 'models/gitm/work']);
% [~,temp] = unix('grep "PBS_JOBNAME=" < pbs_champ_mult_imf_nof_imac.sh | awk ''{print $1}'''); i=find(temp=='='); temp=temp((i+1):(end-1)); %get where this will actually will run
% h=['/nobackup/morozova/' temp '/']; %where did it run? "home" directory
h=[h0 'models/'];
p=[h0 'models/gitm/work/']; %"parent" directory - always come back here after wondering away

cd([h 'gitm/GITM2/src/'])
[~,nLons] = unix('grep "nLons" < ModSize.f90 | awk ''{print $6}'''); nLons=str2double(nLons); % number of longitudes PER block (from ModSize.f90)
[~,nLats] = unix('grep "nLats" < ModSize.f90 | awk ''{print $6}'''); nLats=str2double(nLats); % number of latitudes  PER block
[~,nAlts] = unix('grep "nAlts" < ModSize.f90 | awk ''{print $6}'''); nAlts=str2double(nAlts); % number of altitudes  total (since there is no vertical blocking)

cd([h 'gitm/GITM2/'])
[~,nBlocksLon] = unix('grep "lons" < UAM.in.truncated | awk ''{print $1}'''); nBlocksLon=str2double(nBlocksLon); % number of blocks longitude-wise (from UAM.in)
[~,nBlocksLat] = unix('grep "lats" < UAM.in.truncated | awk ''{print $1}'''); nBlocksLat=str2double(nBlocksLat); % number of blocks latitude -wise
cd(p)

cd(['~/restartOUT']) %what restart files do you want to plot?
nSpecies = 5;       % number of "important" neutral species (from ModEarth.f90)
nSpeciesTotal = 11; % number of    "all"    neutral species
nIons = 10;         % number of    "all"    ion     species

nGhost = 2; %number of ghost cells

%iterate over block files - start in SW and go E first and N second.
for b=1:nBlocksLon*nBlocksLat 
    bb=[num2str(zeros(1,4-length(num2str(b)))) num2str(b)];
    bb(bb==' ')='';
    fn=['b' bb '.rst']
    id=fopen(fn,'r'); %open the current block file
    
    %Lon
    nv=fread(id,1,'int32')/8; %how long is the array to be read?, should be nLons+2*nGhost
    Lon(:,b)=fread(id,nv,'float64'); %n, b - order of indecies in Lon
    nv=fread(id,1,'int32')/8; %just to doublecheck - how many values did you just read?
                              % (this call needs to be here as to keep the cursor position in the file in the correct place - just a FORTRAN convention)    
    %Lat
    nv=fread(id,1,'int32')/8;
    Lat(:,b)=fread(id,nv,'float64'); %t, b
    nv=fread(id,1,'int32')/8;
    
    %Alt
    nv=fread(id,1,'int32')/8;
    Alt(:,b)=fread(id,nv,'float64'); %a, b
    nv=fread(id,1,'int32')/8;
    
    %NeutralDensityS
    for s=1:nSpeciesTotal
        nv=fread(id,1,'int32')/8;
        %how indecies are ordered in NDS : n,t,a,s,b
        NDSt=fread(id,nv,'float64'); % "t" stands for temporary
        NDS(:,:,:,s,b)=reshape(NDSt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost);
        nv=fread(id,1,'int32')/8;
    end
    
    %IonDensityS
    for i=1:nIons
        nv=fread(id,1,'int32')/8;
        %n,t,a,i,b
        IDSt=fread(id,nv,'float64');
        IDS(:,:,:,i,b)=reshape(IDSt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost);
        nv=fread(id,1,'int32')/8;
    end
    
    %Temperature
    nv=fread(id,1,'int32')/8;
    %n,t,a,b
    Tempt=fread(id,nv,'float64');
    Temp(:,:,:,b)=reshape(Tempt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost);
    nv=fread(id,1,'int32')/8;
    
    %IonTemperature
    nv=fread(id,1,'int32')/8;
    %n,t,a,b
    ITempt=fread(id,nv,'float64');
    ITemp(:,:,:,b)=reshape(ITempt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost);
    nv=fread(id,1,'int32')/8;
    
    %electronTemperature
    nv=fread(id,1,'int32')/8;
    %n,t,a,b
    eTempt=fread(id,nv,'float64');
    eTemp(:,:,:,b)=reshape(eTempt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost);
    nv=fread(id,1,'int32')/8;
    
    %Velocity
    nv=fread(id,1,'int32')/8;
    %n,t,a,3,b
    Vt=fread(id,nv,'float64');
    V(:,:,:,:,b)=reshape(Vt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost, 3);
    nv=fread(id,1,'int32')/8;
    
    %IonVelocity
    nv=fread(id,1,'int32')/8;
    %n,t,a,3,b
    IVt=fread(id,nv,'float64');
    IV(:,:,:,:,b)=reshape(IVt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost, 3);
    nv=fread(id,1,'int32')/8;
    
    %VerticalVelocity
    nv=fread(id,1,'int32')/8;
    %n,t,a,s,b
    VVt=fread(id,nv,'float64');
    VV(:,:,:,:,b)=reshape(VVt, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost, nSpecies);
    nv=fread(id,1,'int32')/8 %  Printed to screen just for doublechecking purposes. 
                             %+ It should be an integer, of reasonable size
                             %+ (45630 for stadard 9x9x2x2x50 setup), but anything like 2e10 or -1e-9.7834 is bad. 
    
%F107 - not available in standard GITM distributions! Comment it out.
    nv=fread(id,1,'int32')/8;
    %n,t,a,s,b
    F107t=fread(id,nv,'float64');
    F107(b)=F107t;
    nv=fread(id,1,'int32')/8; %  Printed to screen just for doublechecking purposes. 
                             %+ It should be an integer, of reasonable size
                             %+ (45630 for stadard 9x9x2x2x50 setup), but anything like 2e10 or -1e-9.7834 is bad. 
    
       
%Rho - not available in standard GITM distributions! Comment it out.
    nv=fread(id,1,'int32')/8;
    %n,t,a,s,b
    Rhot=fread(id,nv,'float64');
    Rho(:,:,:,b)=reshape(Rhot, nLons+2*nGhost, nLats+2*nGhost, nAlts+2*nGhost);
    nv=fread(id,1,'int32')/8 %  Printed to screen just for doublechecking purposes. 
                             %+ It should be an integer, of reasonable size
                             %+ (45630 for stadard 9x9x2x2x50 setup), but anything like 2e10 or -1e-9.7834 is bad. 
    
    
    fclose(id); %close the current block file
end

%% 2.PUTTING-TOGETHER section (and removing ghost cells)

% "T" on the end stands for Total or Together
AltT=Alt(2+(1:nAlts),1); %Since there are no blocks in altitude, take it from the first block


% Initialize the variables as big fields of zeros and then paste the values
% from multiple files
NDST=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts, nSpeciesTotal);
IDST=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts, nIons);
TempT=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts);
ITempT=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts);
eTempT=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts);
VT=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts, 3);
IVT=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts, 3);
VVT=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts, nSpecies);
RhoT=zeros(nBlocksLon*nLons, nBlocksLat*nLats, nAlts);

% the following "pasting" is probably easier done with "reshape", but oh well
LatT=[];
for bLt=1:nBlocksLat     % go North (over lats) last  (slower-varying loop)
    LatT=[LatT; Lat(2+(1:nLats),1+(bLt-1)*nBlocksLon)];
    
    LonT=[];
    for bLn=1:nBlocksLon % go East  (over lons) first (faster-varying loop)
        LonT=[LonT; Lon(2+(1:nLons),bLn)];
        NDST((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :, :) = squeeze(NDS(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), :, bLn+(bLt-1)*nBlocksLon));
        IDST((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :, :) = squeeze(IDS(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), :, bLn+(bLt-1)*nBlocksLon));
        TempT((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :) = squeeze(Temp(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), bLn+(bLt-1)*nBlocksLon));
        ITempT((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :) = squeeze(ITemp(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), bLn+(bLt-1)*nBlocksLon));
        eTempT((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :) = squeeze(eTemp(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), bLn+(bLt-1)*nBlocksLon));
        VT((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :, :) = squeeze(V(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), :, bLn+(bLt-1)*nBlocksLon));
        IVT((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :, :) = squeeze(IV(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), :, bLn+(bLt-1)*nBlocksLon));
        VVT((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :, :) = squeeze(VV(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), :, bLn+(bLt-1)*nBlocksLon));
        RhoT((bLn-1)*nLons+(1:nLons), (bLt-1)*nLats+(1:nLats), :) = squeeze(Rho(2+(1:nLons), 2+(1:nLats), 2+(1:nAlts), bLn+(bLt-1)*nBlocksLon));
    end
end

% clear Lon Lat Alt Temp Tempt eTemp eTempt ITemp ITempt NDS NDSt IDS IDSt V Vt IV IVt VV VVt

%% MINIMUMS AND MAXIMUMS

%just so you know what your data looks like!
disp(' ')
disp(' ')
disp('MINIMA AND MAXIMA')
disp(' ')
disp(['LonT ', num2str(min(LonT(:)) ), ' ', num2str(max(LonT(:)) ) ])
disp(['LatT ', num2str(min(LatT(:)) ), ' ', num2str(max(LatT(:)) ) ])
disp(['AltT ', num2str(min(AltT(:))), ' ', num2str(max(AltT(:))) ])
disp(['TempT ', num2str(min(TempT(:))), ' ', num2str(max(TempT(:))) ])
disp(['ITempT ', num2str(min(ITempT(:))), ' ', num2str(max(ITempT(:))) ])
disp(['eTempT ', num2str(min(eTempT(:))), ' ', num2str(max(eTempT(:))) ])
disp('_')
for i=1:nSpeciesTotal
    n=NDST(:,:,:,i);
    disp(['NDST_', num2str(i) ' ' num2str(min(n(:))), ' ', num2str(max(n(:))) ])
end
disp('_')
for i=1:nIons
    n=IDST(:,:,:,i);
    disp(['IDST_', num2str(i) ' ' num2str(min(n(:))), ' ', num2str(max(n(:))) ])
end
disp('_')
disp(['VT ', num2str(min(VT(:))), ' ', num2str(max(VT(:))) ])
disp(['IVT ', num2str(min(IVT(:))), ' ', num2str(max(IVT(:))) ])
disp(['VVT ', num2str(min(VVT(:))), ' ', num2str(max(VVT(:))) ])
disp(['RhoT ', num2str(min(RhoT(:))), ' ', num2str(max(RhoT(:))) ])




%% 3.PLOTTING section

figure(1); %temp on Lat vs Lon
imagesc(LonT*180/pi, LatT*180/pi, squeeze(TempT(:,:,round(nAlts/2))'))
colorbar
set(gca,'YDir','normal')
set(gca,'YTick',LatT*180/pi)
set(gca,'XTick',LonT*180/pi)
title(['Temperature in K at ' num2str(round(AltT(round(nAlts/2))/1000)) ' km.'])
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')


figure(2); %Velocity on Lat vs Lon
contourf(LonT*180/pi, LatT*180/pi, squeeze(VT(:,:,round(nAlts/2),2))')
colorbar
set(gca,'YDir','normal')
set(gca,'YTick',LatT*180/pi)
set(gca,'XTick',LonT*180/pi)
title(['Velocity (m/s) at ' num2str(round(AltT(round(nAlts/2))/1000)) ' km.'])
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')


figure(3); %temp on alt vs Lon
imagesc(LonT*180/pi, AltT, squeeze(TempT(:,1,:))')
colorbar
set(gca,'YTick',AltT)
set(gca,'XTick',LonT*180/pi)
set(gca,'YDir','normal')
title(['Temperature in K at ' num2str(LatT(1)) 'deg Lat'])
xlabel('Longitude (deg)')
ylabel('Altitude (m)')

% 3D cube plot - color is the temperature
figure(4); %temp on alt vs Lat vs Lon
[n,t,a]=meshgrid(LonT*180/pi, LatT*180/pi, AltT/1000);
h=slice(n,t,a,permute(TempT,[2,1,3]),LonT*180/pi,LatT*180/pi,AltT/1000);
set(h,'facealpha',0.1,'edgealpha',0.1)
axis vis3d
colorbar
% set(gca,'ZTick',AltT)
title(['Temperature'])
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
zlabel('Altitude (km)')



% TEMPERATURE PLOT ON A SPHERE
figure(5); clf %temp on alt vs Lat vs Lon

Er=6.378100; %radius of the earth in Mega meters (Mm)
[th,ph,r]=meshgrid(LonT, LatT, AltT/10^6+Er); %create a mesh of lon,lat,rho
x=r.*cos(th).*cos(ph); %convert the mesh locations into cartesian grid
y=r.*sin(th).*cos(ph);
z=r.*sin(ph);

x=reshape(x, nBlocksLat*nLats, nBlocksLon*nLons*nAlts); % since surf doesn't like 3D arrays for x,y and z - squish them
y=reshape(y, nBlocksLat*nLats, nBlocksLon*nLons*nAlts);
z=reshape(z, nBlocksLat*nLats, nBlocksLon*nLons*nAlts);
c=reshape(permute(TempT,[2,1,3]), nBlocksLat*nLats, nBlocksLon*nLons*nAlts); %Temperature has to be permuted because of surf

%the atmosphere plot
h=surf(x,y,z,c);
set(h,'facealpha',1/nAlts,'edgealpha',1/nAlts)

hold on

    load('topo.mat','topo','topomap1');
    % c=(topo<1000) .* (topo>-1000); % a way to get a binary mask showing most
    % % of the continents

    % % One way to plot Earth - just use SPHERE3D written by J De Freitas:
    % % http://www.mathworks.com/matlabcentral/fileexchange/8585-3d-plot-on-a-spherical-surface
    % topo=topo+min(topo(:)); % first, make the topo on the same colorlevel as
    % % temperature (not very smart, but oh well).
    % topo=topo/max(topo(:));
    % topo=topo*(max(c(:))-min(c(:)))+min(c(:));
    % sphere3d(flipud(topo),0,2*pi,-pi/2,pi/2,Er,4,'contour'); % see "help sphere3d"

    C = contourc(topo,[0 0]); % calculate the contour-pairs
    nLevel = 0;               % initialize the counter

    % retrieve contour coordinates (lon,lat) from contourc function output
    while 1==1
        num = C(2,1);        % number of (lon,lat) contour pairs in current contour-piece
        nLevel = nLevel + 1; % index of the current contour-piece
        lon(1:num,nLevel) = C(1,2:num+1)*pi/180; % extract the longitudes in rads
        lat(1:num,nLevel) = C(2,2:num+1)*pi/180-pi/2; % extract the latitudes in rads
        if size(C,2)<num+2; break; end % see if we read all of C
        C = C(:,num+2:end);  % discard what we have already read into lon and lat
    end

    % Replace zeros with NaNs in x and y so that the zeros aren't plotted
    lon(lon == 0) = NaN;
    lat(lat == 0) = NaN;

    %Convert contour coords from polar (lon, lat, rho) into cartesian grid (x,y,z)
    x=Er*cos(lon).*cos(lat);
    y=Er*sin(lon).*cos(lat);
    z=Er*sin(lat);

    %plot the actual contours on the sphere
    l = line(x,y,z);
    set(l,'LineWidth',2,'color',[0 0 1]);

    % % Other ways to draw the Earth
    % [x,y,z]=sphere(50); x=Er*x; y=Er*y; z=Er*z;
    %     h=surf(x,y,z);
    %     set(h,'CData',topo,'FaceColor','texturemap')
    %     set(h,'facealpha',0.9,'edgealpha',0.1')
    % % [x,y,z] = sphere(50);
    % % props.Cdata = topo;
    % % surface(x,y,z,props);

hold off

% view(3)
axis equal %so that the units on all axis have the same visual length
axis vis3d %so the plot-box doesn't change size during rotation
title('Temperature')
xlabel('x (Mm)')
ylabel('y (Mm)')
zlabel('z (Mm)')
colorbar

%% if you have IMAGESCN in your path (code Alexey wrote for 2D-stretched-grid plots)

figure(3); %temp on alt vs Lon
imagescn(LonT*180/pi, AltT, squeeze(TempT(:,1,:))')
colorbar
set(gca,'YTick',AltT)
set(gca,'XTick',LonT*180/pi)
set(gca,'YDir','normal')
title(['Temperature in K at ' num2str(LatT(1)) 'deg Lat'])
xlabel('Longitude (deg)')
ylabel('Altitude (m)')


%% 4. MAT2RST ones

% This code produces new restart files (named cXXXX.rst) and puts a constant
% number in all fields except Lon, Lat, Alt. Useful for diagnosing/troubleshooting
% /writing the code as far as what gets written/read where,
% but NOT(!) for actual initializing GITM or running GITM for real.

q=3.3; %what to put into the restart files: 0 means all values are 0 (except Lon, Lat and Alt), 1.1 means all are 1.1, etc

for b=1:nBlocksLon*nBlocksLat %iterate over block files - start in SW and go E first and N second.
    bb=[num2str(zeros(1,4-length(num2str(b)))) num2str(b)];
    bb(bb==' ')='';
    id=fopen(['c' bb '.rst'],'w');
    
    %Lon
    nv=nLons+2*nGhost;
    fwrite(id,8*nv,'int32');  %how long is the array to be written?, should be nLons+2*nGhost
    fwrite(id,Lon(:,b),'float64'); 
    fwrite(id,8*nv,'int32');  %just to doublecheck - how many values did you just write?
    
    %Lat
    nv=nLats+2*nGhost;
    fwrite(id,8*nv,'int32');  
    fwrite(id,Lat(:,b),'float64'); 
    fwrite(id,8*nv,'int32'); 
    
    %Alt
    nv=nAlts+2*nGhost;
    fwrite(id,8*nv,'int32'); 
    fwrite(id,Alt(:,b),'float64'); 
    fwrite(id,8*nv,'int32'); 
    
    %NDS
    for s=1:nSpeciesTotal
        nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost);
        fwrite(id,8*nv,'int32'); 
        fwrite(id,q*ones(nv,1),'float64'); 
        fwrite(id,8*nv,'int32'); 
    end
    
    %IDS
    for i=1:nIons
        nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost);
        fwrite(id,8*nv,'int32'); 
        fwrite(id,q*ones(nv,1),'float64'); 
        fwrite(id,8*nv,'int32'); 
    end
    
    %Temp
    nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost);
    fwrite(id,8*nv,'int32');
    fwrite(id,q*ones(nv,1),'float64');
    fwrite(id,8*nv,'int32');
    
    %ITemp
    nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost);
    fwrite(id,8*nv,'int32');
    fwrite(id,q*ones(nv,1),'float64');
    fwrite(id,8*nv,'int32');
    
    %eTemp
    nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost);
    fwrite(id,8*nv,'int32');
    fwrite(id,q*ones(nv,1),'float64');
    fwrite(id,8*nv,'int32');
    
    %Velocity
    nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost)*3;
    fwrite(id,8*nv,'int32');
    fwrite(id,q*ones(nv,1),'float64');
    fwrite(id,8*nv,'int32');
    
    %IVelocity
    nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost)*3;
    fwrite(id,8*nv,'int32');
    fwrite(id,q*ones(nv,1),'float64');
    fwrite(id,8*nv,'int32');
    
    %VerticalVelocity
    nv=(nLons+2*nGhost)*(nLats+2*nGhost)*(nAlts+2*nGhost)*nSpecies;
    fwrite(id,8*nv,'int32');
    fwrite(id,q*ones(nv,1),'float64');
    fwrite(id,8*nv,'int32');
    
    
    fclose(id); %close the current block restart file
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
