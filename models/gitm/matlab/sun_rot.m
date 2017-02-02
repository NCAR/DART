%% Simulator of Earth rotating about the Sun, CHAMP rotating about the Earth and showing density above the earth
% this file cane make movies or stills, as in figure 1a in gitm/GITM2/srcDoc/thermo.pdf

% DART $Id$
% CREDIT: Alexey Morozov

%% setup
clc
format compact

f=0; %0 for earth fixed, 1 for star fixed
real_d=1; %(only affects f=1) real distances (1) or shortened (0), just for illustration
noon=~f; %artificially rotate everything so that noon is at the top of the figure?
mov=0; %save a movie?
sun=f; %show sun (distance and size not to scale)?
link=f; %show sun-earth link?
map=1; %show earth map?
champ=1; %show champ?
grace=1; %show grace?
sp=1; %show subsolar point?
aa=1; %show ann arbor?

Re=6.378100; %radius of the earth in Mega meters (Mm)
Rs=2*Re; %sun raidus
if real_d; Rs=695.5/Re*Re; end %the true one, about 109 earth radii

xc=10*Re; %initial location of the earth (center of sun to center of earth)
if real_d; xc=150000/Re*Re; end %the true one (avg), about 24000 earth radii away
yc=0;
zc=0;
rca_c=[xc';yc';zc'];

wba=2*pi/(7*24*60); %rate of sun rotation about its axis in rad/min
wca=2*pi/(365*24*60); %rate of Earth rotation about the sun in rad/min
thdc=(23+26/60)*pi/180; %offset angle of earth roation axis (kd) from the normal vector of the Sun-Earth plane (kc), radians
wed=2*pi/(24*60); %rate of Earth rotation about its axis of rotation in rad/min

%% PRECALCULATE (rotate things into the right frame at every step)
% set missing variables
nLons=36;
nLats=18;
nAlts=50;
nBlocksLon=2;
nBlocksLat=2;

% map precalc
load('topo.mat','topo','topomap1');

C = contourc(topo,[0 0]); % calculate the contour-pairs
nLevel = 0;               % initialize the counter
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

x=Re*cos(lon).*cos(lat);
y=Re*sin(lon).*cos(lat);
z=Re*sin(lat);

rme_e=[x(:)';y(:)';z(:)'];

% atmosphere
[th,ph]=meshgrid(LonT*pi/180, LatT*pi/180); %create a mesh of lon,lat
r=Re;%+AltT(walt)/10^6;
x=r*cos(th).*cos(ph); %convert the mesh locations into cartesian grid
y=r*sin(th).*cos(ph);
z=r*sin(ph);

xe=x(:);
ye=y(:);
ze=z(:);
rte_e=[xe';ye';ze'];

%% rotation precalc

% % only plot champ data when truth data is available
% i = find( tc>0 & tc<nom );
% tc=tc(i);
% LonC=LonC(i);
% LatC=LatC(i);
% AltC=AltC(i);
% RhoC=RhoC(i);

%preallocate
rfe_ds=zeros(3,length(find(tc>0 & tc<1440))); %r is for position vector, fe is postion of "f" with respect to "e" (or vector from point "e" to "f"), "_e" is resolved in frame e, s on the end is for "store"
rge_ds=rfe_ds;
% rme_ds=zeros(3,size(rme_e,2),length(tc)); %450e6, so just rotate map at every step. Must be a better way to store the map
rse_ds=rfe_ds;
rre_ds=rfe_ds;

rfa_as=rfe_ds;
rga_as=rfe_ds;
rsa_as=rfe_ds;
rra_as=rfe_ds;

%grace
LonGC=interp1np(tg,LonG,tc);
LatGC=interp1(tg,LatG,tc);
AltGC=interp1(tg,AltG,tc);
RhoGC=interp1(tg,RhoG,tc);

for i=find(tc>0,1):find(tc<1440,1,'last')

    %champ
    xe=(AltC(i)/10^6+Re).*cos(LonC(i)*pi/180).*cos(LatC(i)*pi/180);
    ye=(AltC(i)/10^6+Re).*sin(LonC(i)*pi/180).*cos(LatC(i)*pi/180);
    ze=(AltC(i)/10^6+Re).*sin(LatC(i)*pi/180);
    rfe_e=[xe';ye';ze'];
    rfe_d=om(3,wed*tc(i)-noon*pi/2)'*rfe_e; rfe_ds(:,i)=rfe_d; %earth rotation about axis
    rfe_c=om(2,thdc)'*rfe_d; %rotation axis tilt
    rfe_a=om(3,wca*(tc(i)+d1*1440))'*rfe_c; %earth rotation about sun
    red_a=[0;0;0]; rdc_a=[0;0;0]; %origins of c, d and e frames are colocated
    rca_a=om(3,wca*(tc(i)+d1*1440))'*rca_c; %earth rotation about sun, transpose to get Oac=Oca'
    rea_a=red_a+rdc_a+rca_a;
    rfa_a=bsxfun(@plus,rfe_a,rea_a); rfa_as(:,i)=rfa_a; %position wrt the sun
    
    %grace
    xe=(AltGC(i)/10^6+Re).*cos(LonGC(i)*pi/180).*cos(LatGC(i)*pi/180);
    ye=(AltGC(i)/10^6+Re).*sin(LonGC(i)*pi/180).*cos(LatGC(i)*pi/180);
    ze=(AltGC(i)/10^6+Re).*sin(LatGC(i)*pi/180);
    rge_e=[xe';ye';ze'];
    rge_d=om(3,wed*tc(i)-noon*pi/2)'*rge_e; rge_ds(:,i)=rge_d; %earth rotation about axis
    rge_c=om(2,thdc)'*rge_d; %rotation axis tilt
    rge_a=om(3,wca*(tc(i)+d1*1440))'*rge_c; %earth rotation about sun
    red_a=[0;0;0]; rdc_a=[0;0;0]; %origins of c, d and e frames are colocated
    rca_a=om(3,wca*(tc(i)+d1*1440))'*rca_c; %earth rotation about sun, transpose to get Oac=Oca'
    rea_a=red_a+rdc_a+rca_a;
    rga_a=bsxfun(@plus,rge_a,rea_a); rga_as(:,i)=rga_a; %position wrt the sun
    
    %map
    %no way - too gigantic of a matrix to precompute, instead done at every step below
    
    %subsolar = sp
    LonSi=mod( (24*60-mod(tc(i),24*60))/(24*60)*360+180 , 360); %subsolar point (point on Earth closest to the Sun) moves Westward at 360deg/day=15deg/hr. It is where the local noon is (ie where it is 12:00 military time, not 00:00).
    LatSi=  (23+26/60)*sin(2*pi/365*((d1+tc(i)/(24*60))-78) ); %subsolar point moves between the Tropic of Cancer (Northern Solstice, let's say Jun 22) and the Tropic of Capricorn (Southern Solstice, let's say Dec 21).
    %+ so it is a sinusoid with magnitude of +-23.43deg (that's where the Tropics are) and frequency of 2pi rev/365 days. The phase shift is approximately datenum(2002,3,20)-datenum(2002,1,1)
    AltSi=393983.5; %about 400km, one of the gridpoints in GITM in this particular setup
    xe=(AltSi/10^6+Re).*cos(LonSi*pi/180).*cos(LatSi*pi/180);
    ye=(AltSi/10^6+Re).*sin(LonSi*pi/180).*cos(LatSi*pi/180);
    ze=(AltSi/10^6+Re).*sin(LatSi*pi/180);
    rse_e=[xe';ye';ze']; nts_es(:,i)=[LonSi;LatSi;AltSi];
    rse_d=om(3,wed*tc(i)-noon*pi/2)'*rse_e; rse_ds(:,i)=rse_d; %earth rotation about axis
    rse_c=om(2,thdc)'*rse_d; %rotation axis tilt
    rse_a=om(3,wca*(tc(i)+d1*1440))'*rse_c; %earth rotation about sun
    red_a=[0;0;0]; rdc_a=[0;0;0]; %origins of c, d and e frames are colocated
    rca_a=om(3,wca*(tc(i)+d1*1440))'*rca_c; %earth rotation about sun, transpose to get Oac=Oca'
    rea_a=red_a+rdc_a+rca_a;
    rsa_a=bsxfun(@plus,rse_a,rea_a); rsa_as(:,i)=rsa_a; %position wrt the sun
    
    %ann arbor = aa
    LonAi=277.5;
    LatAi=42.5;
    AltAi=393983.5; %about 400km, one of the gridpoints in GITM in this particular setup
    xe=(AltAi/10^6+Re).*cos(LonAi*pi/180).*cos(LatAi*pi/180);
    ye=(AltAi/10^6+Re).*sin(LonAi*pi/180).*cos(LatAi*pi/180);
    ze=(AltAi/10^6+Re).*sin(LatAi*pi/180);
    rre_e=[xe';ye';ze'];
    rre_d=om(3,wed*tc(i)-noon*pi/2)'*rre_e; rre_ds(:,i)=rre_d; %earth rotation about axis
    rre_c=om(2,thdc)'*rre_d; %rotation axis tilt
    rre_a=om(3,wca*(tc(i)+d1*1440))'*rre_c; %earth rotation about sun
    red_a=[0;0;0]; rdc_a=[0;0;0]; %origins of c, d and e frames are colocated
    rca_a=om(3,wca*(tc(i)+d1*1440))'*rca_c; %earth rotation about sun, transpose to get Oac=Oca'
    rea_a=red_a+rdc_a+rca_a;
    rra_a=bsxfun(@plus,rre_a,rea_a); rra_as(:,i)=rra_a; %position wrt the sun
    
end




%% PLOT
close
ift=1; %is first time?

for i=find(tc>0,1)+195; %50:290%1:length(tc)%160:161 (1:240)+80 50:290 %195
    figure(1)
    clf
%     subplot(1,2,1)
%     subplota(1,2,1,1,.02,0)
% axes('position', [.07 .1 .4 .8] );

    hold on
    
    %     disp(round(i/length(tc)*100))
    disp(i)
    %% SUN rotation simulator
    if sun
        [th,ph]=meshgrid(LonT*pi/180, LatT*pi/180); %create a mesh of lon,lat
        r=Rs;
        x=r*cos(th).*cos(ph); %convert the mesh locations into cartesian grid
        y=r*sin(th).*cos(ph);
        z=r*sin(ph);
        
        xb=x(:);
        yb=y(:);
        zb=z(:);
        rb=[xb';yb';zb'];
        ra=om(3,-wba*tc(i))*rb;
        
        x=reshape(ra(1,:), nBlocksLat*nLats, nBlocksLon*nLons); % since surf doesn't like 3D arrays for x,y and z - squish them
        y=reshape(ra(2,:), nBlocksLat*nLats, nBlocksLon*nLons);
        z=reshape(ra(3,:), nBlocksLat*nLats, nBlocksLon*nLons);
        
        %     c=reshape(permute(log10(RhoT(:,:,35,j)) ,[2,1,3]), nBlocksLat*nLats, nBlocksLon*nLons*1); %Temperature has to be permuted because of surf
        c=permute(log10(RhoT(:,:,35,1)) ,[2,1,3]); %Temperature has to be permuted because of surf
        x=[x x(:,1)];
        y=[y y(:,1)];
        z=[z z(:,1)];
        c=-11.42+0*[c c(:,1)];
        h=surf(x,y,z,c);
%         set(h,'facealpha',1/2,'edgealpha',0)%,'FaceLighting','phong','FaceColor',[1 2/3 0],'AmbientStrength',0.9);light('Position',[1 0 0],'Style','infinite');
set(h,'edgecolor','none')%,'FaceLighting','phong','FaceColor',[1 2/3 0],'AmbientStrength',0.9);light('Position',[1 0 0],'Style','infinite');
    end
    
    %% Sun-earth link
    if link
        rca_a=om(3,wca*(tc(i)+d1*1440))'*rca_c;
        line([0 rca_a(1)],[0 rca_a(2)],[0 rca_a(3)],'color','k')
    end
    %% map
    if map
        
        %Convert contour coords from polar (lon, lat, rho) into cartesian grid (x,y,z)
        rme_d=om(3,wed*tc(i)-noon*pi/2)'*rme_e; %earth rotation about axis
        rme_c=om(2,thdc)'*rme_d; %rme_ds(:,:,i)=rme_c; %rotation axis tilt
        rme_a=om(3,wca*(tc(i)+d1*1440))'*rme_c; %earth rotation about sun
        red_a=[0;0;0]; rdc_a=[0;0;0]; %origins of c, d and e frames are colocated
        rca_a=om(3,wca*(tc(i)+d1*1440))'*rca_c; %earth rotation about sun, transpose to get Oac=Oca'
        rea_a=red_a+rdc_a+rca_a;
        rma_a=bsxfun(@plus,rme_a,rea_a);
        
        
        %plot the actual contours on the sphere
        if f==0
            l = line(rme_d(1,:),rme_d(2,:),rme_d(3,:));
        else
            l = line(rma_a(1,:),rma_a(2,:),rma_a(3,:));
        end
        set(l,'LineWidth',1,'color','k');
        
    end
    
    %% champ
    if champ
        if f==0
            ra=rfe_ds;
        else
            ra=rfa_as;
        end
        %plot the full traj
        l = line(ra(1,:),ra(2,:),ra(3,:));
        set(l,'LineWidth',1,'color','k');
        
        %plot the current dot
        h_f = plot3(ra(1,i),ra(2,i),ra(3,i),...
            'ro','LineWidth',2,'MarkerFaceColor',[255 153 153]/255,...
            'markersize',10);
    end
    
    %% grace
    if grace
        if f==0
            ra=rge_ds;
        else
            ra=rga_as;
        end
        %plot the full traj
        l = line(ra(1,:),ra(2,:),ra(3,:));
        set(l,'LineWidth',1,'color','k');
        
        %plot the current dot
        h_g = plot3(ra(1,i),ra(2,i),ra(3,i),...
            's','color',[1 .5 0],...
            'LineWidth',2,'MarkerFaceColor','y',...
            'markersize',10);
    end
    
    %% sp
    if sp
        if f==0
            ra=rse_ds;
        else
            ra=rsa_as;
        end
        %plot the full traj
%         l = line(ra(1,:),ra(2,:),ra(3,:));
        set(l,'LineWidth',1,'color','k');
        %plot the current dot
        h_s = plot3(ra(1,i),ra(2,i),ra(3,i),...
            'b^','LineWidth',2,'MarkerFaceColor',[135 206 250]/255,...
            'markersize',10);
    end
    
    %% aa
    if aa
        if f==0
            ra=rre_ds;
        else
            ra=rra_as;
        end
        %plot the full traj
%         l = line(ra(1,:),ra(2,:),ra(3,:));set(l,'color','k');
        %plot the current dot
        h_r = plot3(ra(1,i),ra(2,i),ra(3,i),...
            'gd','color',[0 0.5 0],'LineWidth',2,'MarkerFaceColor',[0 1 0],...
            'markersize',10);
        
        
    end
    
    
    %% atmosphere
    
    rte_d=om(3,wed*tc(i)-noon*pi/2)'*rte_e; %earth rotation about axis
    rte_c=om(2,thdc)'*rte_d; %rotation axis tilt
    rte_a=om(3,wca*(tc(i)+d1*1440))'*rte_c; %earth rotation about sun
    red_a=[0;0;0]; rdc_a=[0;0;0]; %origins of c, d and e frames are colocated
    rca_a=om(3,wca*(tc(i)+d1*1440))'*rca_c; %earth rotation about sun, transpose to get Oac=Oca'
    rea_a=red_a+rdc_a+rca_a;
    rta_a=bsxfun(@plus,rte_a,rea_a);
    
    if f==0
        ra=rte_d;
    else
        ra=rta_a;
    end
    
    x=reshape(ra(1,:), nBlocksLat*nLats, nBlocksLon*nLons); % since surf doesn't like 3D arrays for x,y and z - squish them
    y=reshape(ra(2,:), nBlocksLat*nLats, nBlocksLon*nLons);
    z=reshape(ra(3,:), nBlocksLat*nLats, nBlocksLon*nLons);
    
    j=find(tt<tc(i),1,'last');
    if isempty(j)
        j=1;
    end
    %     c=reshape(permute(log10(RhoT(:,:,35,j)) ,[2,1,3]), nBlocksLat*nLats, nBlocksLon*nLons*1); %Temperature has to be permuted because of surf
    c=permute(10^12*RhoT(:,:,35,j) ,[2,1,3]); %Temperature has to be permuted because of surf
    x=[x x(:,1)];
    y=[y y(:,1)];
    z=[z z(:,1)];
    c=[c c(:,1)];
    h_atm=surf(x,y,z,c);
    set(h_atm,'facealpha',1/2,'edgealpha',0)
    
    
    %% misc
    %     grid on
    %     box on
    %     set(gca,'position',[.1 .06 .7 0.9])
    %     view(3)
    %axis equal %so that the units on all axis have the same visual length
    %     axis vis3d %so the plot-box doesn't change size during rotation
    
    [yy mo dd hh mm ss]=datevec(datenum(2000+y2,0,0)+d1+tc(i)/1440);
    %     title(['Mass density (kg/m3) at 400 km, ' ...
    %         num2str(yy) '-' num2str(mo,'%02.0f') '-' num2str(dd,'%02.0f')...
    %         ' ' num2str(hh,'%02.0f') ':' num2str(mm,'%02.0f') ':' num2str(ss,'%02.0f')...
    %         ' UTC.'])
    
    % text(Re,-Re, -Re, [num2str(hh,'%02.0f') ':' num2str(mm,'%02.0f') 'UT'])
    
    
    %     xlabel('Longitudonal distance in Sun-Earth frame (Mm)')
    %     ylabel('Transverse distance (Mm)')
    
    
    if f==0
        ylim([-Re+1 Re+1])
        axis off
        th=0:.01:2*pi;
        r=Re+1;
%         h=patch(r*cos(th),r*sin(th),0*th,[1 1 1]);
%         set(h,'facealpha',0,'edgealpha',1)
        l=line(r*cos(th),r*sin(th),0*th);
        set(l,'color','k');
        
        for k=0:23
            th2=k*2*pi/24;
            r2=Re+.7;
            l=line([r r2]*cos(th2),[r r2]*sin(th2),[0 0]);
            set(l,'color','k');
        end
        for k=0:3
            th2=k*2*pi/4;
            r2=Re+1.6;
            text(r2*cos(th2)-.06*Re,r2*sin(th2),0,num2str( mod((k+1)*6,24) ,'%02.0f'),'fontsize',12);
        end
        text(.3*Re,-r2*1.1,0,'(a)','fontsize',12);
        cb=colorbar;
        ylabel(cb,'Mass density  x10^{-12} [kg m^{-3}]','fontsize',12)
        set(cb,'fontsize',12);
        %     y=get(cb,'Ytick');
        %     set(cb,'Ytick',min(y):.1:max(y));
        
        %      set(gca, 'position',[0 0 1*.613/.815 1])
        
         set(gcf,'position',[521   242   700   273]); %for double figure
                                                      %            set(gcf,'position',[521   242   420   273]); %for single figure
        %     set(gcf,'position',[521   242   734   462]); %bigger size
        %     drawnow
        
        l=legend([h_f h_g h_s h_r],...
            {'CHAMP','GRACE','Subsolar pt.','Ann Arbor'},...
            'location','SE');
        
    else
        if ~real_d; xlim([-80 80]); end
        set(gcf,'position',[345 49 911 681]); 
    end
    
%     cl=10^12*10.^[-11.8 -11.2];
cl=[2 7];
    set(gca,'clim',cl)
    set(gcf,'color',[1 1 1])
    axis equal
    hold off
    
    
    
    if mov
        if  ift
            ift=0;
            fr=getframe(gcf);
            [im1] = frame2im(fr);
            [im2,ma]=rgb2ind(im1,128);
            %         imwrite(im2,ma,fname,'LoopCount',0,'DelayTime',1);
            
%             aviobj = avifile('mc.avi','colormap',ma,'compression','msvc','quality',75);
%             aviobj = addframe(aviobj,im2);
            
            
        else
            
            
            fr=getframe(gcf);
            im1 = frame2im(fr);
            im2=rgb2ind(im1,ma);
            %         imwrite(im2,ma,fname,'WriteMode','append');
            aviobj = addframe(aviobj,im2);
            
        end
        
    end
    
end

% aviobj = close(aviobj);
%  set(h_atm,'facealpha',1,'edgealpha',1,'edgecolor','none')
%  
% set(gcf,'renderer','painters')

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
