%% common plot parameters CR

% DART $Id$
% CREDIT: Alexey Morozov


al=0.5; %transparency
scale=10^12; %multiplication factor for densities (I'm plotting RhoT*scale as RhoT was causing matlab plotting problems)
mavg=1.5; %how many hours to average data over
ex='6'; %I had multiple examples and this was #6, affects the name of the generated plot file
pri=1; %save plots into pdfs (1) or just display on screen (0)?
cfh_h=cfh; %what localization to use (cfh comes from pbs_file.sh, but you can plug in your own value for experimentation - 0.6 great circle rad or 0.7 rad or 1.2 or 0.1 etc)

% lon1=20;
% lat1=1;
% alt1=37;
lon1=wlon; %longitude cell index of interest (for plots - I picked 1 location (Ann Arbor, MI)
lat1=wlat; %latitude cell index of interest
alt1=walt; %altitude cell index of interest
alt2=min(alta):max(alta); %altitude cell indeces that the satellite visited
% alt2=30:40;

% LonD(lon1)
% LatD(lat1)
% AltD(alt1)/1000

time_p=[min(td) max(td)];

yl=[0 12]*1e-12*scale;
% yl=[0 5]*1e-12*scale;
tsub=9+[0 3]; psub= [0.7 0.7 .3 .3];


%% Rho sawtooth along subsolar point 
figure(1);clf
set(gcf,'position',[50 50 500 350])


% [t_int,~,i_ts] = intersect(round(tt*60*60),round(ts*60*60) ); %round to seconds 
[tts,rs,gs]=loc(tt, LonTS, LatTS, AltTS, ...
    tt, LonCT, LatCT, AltCT, RhoCT, ...
    cfh_h, cfv);

time_plot([],RhoTS,'GITM simulation', ...
    tm,RhoMS,'GITM without EAKF', ...
    [],RhoTS,'GITM truth simulation', ...
    [],kind_sd+0*tt,'GITM truth simulation +/-SD', ...
    [],rs,50,[ones(size(gs)) 1-gs 1-gs ],'CHAMP measurement localized', ...
    td,RhoDrmS,RhoDomS,'EAKF ensemble mean', ...
    td,RhoDrsS,RhoDosS,'EAKF ensemble mean +/- SD', ...
    [],[],[], ...
    [],[4 1],-1.5:.01:3,'na', ...
    'SE','none',[],al,scale)


ylim(yl)
xlim(time_p)
ylabel('Mass density  x10^{-12} [kg m^{-3}]')
xlabel('(a)  Hours since 00UT 12/1/2002')

set(gca,'XTick',0:3:48)

if pri; set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 .125 8.1 5.7]);  print(gcf,'-dpdf',['c' ex '_sp']); end

%% Rho sawtooth at aa

figure(2);clf
set(gcf,'position',[100 100 500 350])

[tts,rs,gs]=loc(tt, LonTA, LatTA, AltTA, ...
    tt, LonCT, LatCT, AltCT, RhoCT, ...
    cfh, cfv);

time_plot([],RhoTA,'GITM truth simulation', ...
    tm,RhoMA,'GITM without EAKF', ...
    [],RhoTA,'GITM truth simulation', ...
    [],kind_sd,'GITM truth simulation +/- SD', ...
    tts,rs,50,[ones(size(gs)) 1-gs 1-gs ],'CHAMP measurement localized', ...
    td,squeeze( RhoDr(lon1,lat1,alt1,1,:) ),squeeze( RhoDo(lon1,lat1,alt1,1,:) ),'EAKF ensemble mean', ...
    td,squeeze( RhoDr(lon1,lat1,alt1,2,:) ),squeeze( RhoDo(lon1,lat1,alt1,2,:) ),'EAKF ensemble mean +/- SD', ...
    [],[],[], ...
    [],[4 1],-1.5:.01:3,'na', ...
    'NE','none',[],al,scale)%[4 1 2 3]

ylim(yl)
xlim(time_p)
ylabel('Mass density  x10^{-12} [kg m^{-3}]')
set(gca,'XTick',0:3:48)
xlabel('(a)  Hours since 00UT 12/1/2002')
% set(gca,'XTick',0:0.5:2)
% xlabel('Days since 2002-12-01')

if pri; set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 .125 8.1 5.7]);  print(gcf,'-dpdf',['c' ex '_1l']); end



%% Rho sawtooth along champ path 
figure(3);clf
set(gcf,'position',[150 150 500 350])

[tts,rs,gs]=loc(tt, LonTC, LatTC, AltTC, ...
    tt, LonCT, LatCT, AltCT, RhoCT, ...
    cfh_h, cfv);

[t90m, r90m]=m90(tm,RhoMC, mavg);
[t90c, r90c]=m90(tc,RhoC, mavg);
[t90cu, r90cu]=m90(tc,RhoCu, mavg);
[t90d, r90dr]=m90(td,RhoDrmC, mavg);
[t90d, r90do]=m90(td,RhoDomC, mavg);
[t90du, r90dru]=m90(td,RhoDrsC, mavg);
[t90du, r90dou]=m90(td,RhoDosC, mavg);

time_plot([],RhoTC,'GITM truth simulation', ...
    t90m,r90m,'GITM without EAKF', ...
    t90c,r90c,'CHAMP measurement', ...
    t90cu,r90cu,'CHAMP measurement +/- SD', ...
    [],rs,50,[ones(size(gs)) 1-gs 1-gs ],'CHAMP measurement localized', ...
    t90d,r90dr,r90do,'EAKF ensemble mean', ...
    t90du,r90dru,r90dou,'EAKF ensemble mean +/- SD', ...
    [],[],[], ...
    [],[4 1],-1.5:.01:3,'na', ...
    'NE','none',[1 2 5 3 4],al,scale)

% ylim(yl)
ylim([0 8]*1e-12*scale)
xlim(time_p)
ylabel('Mass density  x10^{-12} [kg m^{-3}]')
xlabel('(a)  Hours since 00UT 12/1/2002')

set(gca,'XTick',0:3:48)

% subview(2,[24 27], [1.5 3.5]*1e-12, [.6 .7 .3 .25])
% subview(2,9+[0 3], [0.5 7]*1e-12*scale, [.6 .7 .3 .25]);
% ylim([1 5.5]*1e-12)
% subview(2,[21 24], [2 4]*1e-12, [.6 .7 .3 .25])

if pri; set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 .125 8.1 5.7]);  print(gcf,'-dpdf',['c' ex '_cl']); end

%% Rho RMSPE along champ path 


% i = plot_RMSPE_xxx_evolution_ALEX('obs_diag_output.nc','spread');
% j=i{5}; %height (1 - 0km, 2 - 100km, ..., 5 - 400km)
% te=j(:,1);%1 is time, 2 is RMSPE, 3 is spread, 4 is number of observations possible, 5 is number observations used 
% te=te - datenum([2002 12 01 00 00 00]);   %td=Time EAKF (number of minutes since beginning of 12/01/2002)
% te=te*24; %in hrs

i_td=find(td>=24 & td<48);%calculate RMSPE only over the last day
[t_int,i_tt,~] = intersect(round(tt*60*60),round(td(i_td)*60*60) ); %round to seconds because reading from Netcdf files introduces small roundoff noise. Last output is omitted assuming that td is a subset of tt and also because that output would be indexing the wrong array (clipped td instead of full td)
t_int=t_int/3600; %just if you want to know what the intersection is

rms_ct=sqrt(  mean( (RhoCT-RhoTC ).^2 )  )
rmsn_ct=rms_ct/sqrt(  mean( RhoCT.^2 )  )
rms_c=sqrt(  mean( (RhoCT(i_tt)-RhoDomC(i_td) ).^2 )  )
rmsn_c=rms_c/sqrt(  mean( RhoCT(i_tt).^2 )  )

rms_gt=sqrt(  mean( (RhoGT-RhoTG ).^2 )  ) %
rmsn_gt=rms_gt/sqrt(  mean( RhoGT.^2 )  )
rms_gd=sqrt(  mean( (RhoGT(i_tt)-RhoDrmG(i_td) ).^2 )  ) %
rmsn_gd=rms_gt/sqrt(  mean( RhoGT(i_tt).^2 )  )




[t90tt, r90tt]=m90(tt, 100*abs(RhoTC-RhoCT)./RhoCT, mavg);

[t90t, r90t]=m90(tm, 100*abs(RhoMC-RhoCT)./RhoCT, mavg);

[~,i_tt,i_td] = intersect(round(tt*60*60),round(td*60*60) ); %round to seconds because reading from Netcdf files introduces small roundoff noise. Last output is omitted assuming that td is a subset of tt and also because that output would be indexing the wrong array (clipped td instead of full td)
[t90d, r90d]=m90(td(i_td), 100*abs(RhoCT(i_tt)-RhoDomC(i_td))./RhoCT(i_tt), mavg);



figure(31);clf
set(gcf,'position',[150 150 500 350])

time_plot([],r90tt,'GITM truth simulation', ...
    t90t,r90t,'GITM without EAKF', ...
    [],RhoC,'CHAMP measurement', ...
    [],RhoC,'CHAMP measurement', ...
    [],1,50,[],'na', ...
    t90d,r90d,r90d,'EAKF posterior', ...
    [],[],[],'EAKF ensemble mean +/- SD', ...
    [],[],[],...te,100*j(:,5)./j(:,4),'Percent of observations used', ...
    [],[4 1],-1.5:.01:3,'na', ...
    'NE','none',[],al,1)

ylim([0 100])
xlim(time_p)
ylabel('Absolute percentage error along CHAMP path')
xlabel('(c)  Hours since 00UT 12/1/2002')

set(gca,'XTick',0:3:48)

text(18,45,['RMSPE (2nd day) along CHAMP path = ' num2str(round(100*rmsn_c) ) '%'],'FontSize', 12 )

% subview(2,[21 24], [1.5 6.75]*1e-12, [.6 .7 .3 .25])


% ylim([1 5.5]*1e-12)
% subview(2,[21 24], [2 4]*1e-12, [.6 .7 .3 .25])

if pri; set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 .125 8.1 5.7]);  print(gcf,'-dpdf',['c' ex '_cr']); end




%% Rho sawtooth along grace path 
figure(32);clf
set(gcf,'position',[150 150 500 350])

[tts,rs,gs]=loc(tt, LonTG, LatTG, AltTG, ...
    tt, LonCT, LatCT, AltCT, RhoCT, ...
    cfh_h, cfv);

[t90m, r90m]=m90(tm,RhoMG, mavg);
[t90g, r90g]=m90(tg,RhoG, mavg);
[t90gu, r90gu]=m90(tc,RhoGu, mavg);
[t90d, r90dr]=m90(td,RhoDrmG, mavg);
[t90d, r90do]=m90(td,RhoDomG, mavg);
[t90du, r90dru]=m90(td,RhoDrsG, mavg);
[t90du, r90dou]=m90(td,RhoDosG, mavg);

time_plot([],RhoTG,'GITM truth simulation', ...
    t90m,r90m,'GITM without EAKF', ...
    t90g,r90g,'GRACE measurement', ...
    t90gu,r90gu,'GRACE measurement +/- SD', ...
    [],rs,50,[ones(size(gs)) 1-gs 1-gs ],'CHAMP measurement localized', ...
    t90d,r90dr,r90do,'EAKF ensemble mean', ...
    t90du,r90dru,r90dou,'EAKF ensemble mean +/- SD', ...
    [],[],[], ...
    [],[4 1],-1.5:.01:3,'na', ...
    'NE','none',[1 2 5 3 4],al,scale)

ylim([0 3]*1e-12*scale)
xlim(time_p)
ylabel('Mass density  x10^{-12} [kg m^{-3}]')
xlabel('(b)  Hours since 00UT 12/1/2002')

set(gca,'XTick',0:3:48)

% subview(2,[24 27], [1.5 3.5]*1e-12, [.6 .7 .3 .25])
% subview(2,tsub, [0 4]*1e-12*scale, [.6 .5 .3 .25]);
% ylim([1 5.5]*1e-12)
% subview(2,[21 24], [2 4]*1e-12, [.6 .7 .3 .25])

if pri; set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 .125 8.1 5.7]);  print(gcf,'-dpdf',['c' ex '_gl']); end


%% Rho RMSPE along grace path 

i_td=find(td>=24 & td<48);%calculate RMSPE only over the last day
[t_int,i_tt,~] = intersect(round(tt*60*60),round(td(i_td)*60*60) ); %round to seconds because reading from Netcdf files introduces small roundoff noise. Last output is omitted assuming that td is a subset of tt and also because that output would be indexing the wrong array (clipped td instead of full td)
t_int=t_int/3600; %just if you want to know what the intersection is

rms_g=sqrt(  mean( (RhoGT(i_tt)-RhoDomG(i_td) ).^2 )  );
rmsn_g=rms_g/sqrt(  mean( RhoGT(i_tt).^2 )  )

[t90tt, r90tt]=m90(tt, 100*abs(RhoTG-RhoGT)./RhoGT, mavg);

[t90t, r90t]=m90(tm, 100*abs(RhoMG-RhoGT)./RhoGT, mavg);

[~,i_tt,i_td] = intersect(round(tt*60*60),round(td*60*60) ); %round to seconds because reading from Netcdf files introduces small roundoff noise. Last output is omitted assuming that td is a subset of tt and also because that output would be indexing the wrong array (clipped td instead of full td)
[t90d, r90d]=m90(td(i_td), 100*abs(RhoGT(i_tt)-RhoDomG(i_td))./RhoGT(i_tt), mavg);



figure(33);clf
set(gcf,'position',[150 150 500 350])

time_plot([],r90tt,'GITM truth simulation', ...
    t90t,r90t,'GITM without EAKF', ...
    [],RhoC,'CHAMP measurement', ...
    [],RhoC,'CHAMP measurement', ...
    [],1,50,[],'na', ...
    t90d,r90d,r90d,'EAKF posterior', ...
    [],[],[],'EAKF ensemble mean +/- SD', ...
    [],[],'Percent of observations used', ...
    [],[4 1],-1.5:.01:3,'na', ...
    'NE','none',[],al,1)

ylim([0 100])
xlim(time_p)
ylabel('Absolute percentage error along GRACE path')
xlabel('(d)  Hours since 00UT 12/1/2002')

set(gca,'XTick',0:3:48)

text(16,25,['RMSPE (2nd day) along GRACE path = ' num2str(round(100*rmsn_g) ) '%'],'FontSize', 12 )

% subview(2,[21 24], [1.5 6.75]*1e-12, [.6 .7 .3 .25])


% ylim([1 5.5]*1e-12)
% subview(2,[21 24], [2 4]*1e-12, [.6 .7 .3 .25])

if pri; set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 .125 8.1 5.7]);  print(gcf,'-dpdf',['c' ex '_gr']); end


%% F10.7 sawtooth NOCHAMP

figure(5);clf
set(gcf,'position',[250 250 500 350])

time_plot(tf1,(f107T+f107Ta)/2,'F10.7 as measured by NOAA', ...
    tm,130+0*tm,'F10.7 used in GITM without EAKF', ...
    [],RhoC,'na', ...
    [],[],'na', ...
    [],1,50,[],'na', ...
    td,squeeze( f107Dr(1,1,:) ),squeeze( f107Do(1,1,:) ),'EAKF ensemble mean', ...
    td,squeeze( f107Dr(1,2,:) ),squeeze( f107Do(1,2,:) ),'EAKF ensemble mean +/- SD', ...
    [],[],[], ...
    [130 25],[3 0],100:1:230,'Initial F10.7 ensemble distribution', ...
    'NE','none',[],al,1)

% ylim([-1 3])
xlim(time_p)
ylabel('F_{10.7} [SFU]')
set(gca,'XTick',0:3:48)
xlabel('Hours since 00UT 12/1/2002')
% set(gca,'XTick',0:0.5:2)
% xlabel('Days since 2002-12-01')

if pri; set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 .125 8.1 5.7]);  print(gcf,'-dpdf',['c' ex '_f']); end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
