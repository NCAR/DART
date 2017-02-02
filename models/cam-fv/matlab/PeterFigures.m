% PeterFigures generates 'publication-quality' figures for Peter Lauritzen
%
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% the figures Peter has in his noise paper are:
%
% bluefire:/ptmp/raeder/Single_CAM/Halo-diff_12-5-2000/...
% undamped.VS.266hPa.ps         from Prior_memb1_VS_266.nc
% damped.VS.266hPa.ps           from caminput_VS_266.nc
% undamped-damped.VS.266hPa.ps  from 1098_Prior-caminput.nc
%
% Please confirm that they look like what's in the paper Jeff has.
%
% Prior_memb1_VS_266.nc  == /ptmp/raeder/POP_force/POP12/obs_1098/Prior_Diag.nc
%                             time 1, member 1(copy 3), VS, level 266 hPa(k=14)
% caminput_VS_266.nc     == Halo-diff_12-5-2000/caminput.nc  only 1 timeslot there.
% 1098_Prior-caminput.nc == ncdiff Prior_memb1_VS_266.nc caminput_VS_266.nc
%
% Peter would like to have noise at mid-latitudes highlighted, 
% as well as the polar stuff.
% I've higlighted some candidates in my talk:
% mirage3:/ptmp/raeder/friday_talk_init_eval_IPCC.pptx

%% Get the colormap that matches ncview ...

gauss3 = load('/fs/image/home/thoar/matlab/3gauss.ncmap') ./ 255;
bob    = gauss3(1:1:256,:);

%% Get all the data so we can determine consistent limits
% TJH - all the other figures in the paper are [-180,180], which
% is not what we have. Should swap halves to be consistent,
% and repeat the wrap longitude so we can label [-180,180]

% Undamped ... the noisy one.

undamped.fname = '/gpfs/ptmp/raeder/Single_CAM/Halo-diff_12-5-2000/Prior_memb1_VS_266.nc';
undamped.lon   = nc_varget(undamped.fname,'slon');
undamped.lat   = nc_varget(undamped.fname,'lat');
undamped.lev   = nc_varget(undamped.fname,'lev');
undamped.time  = nc_varget(undamped.fname,'time');
myVS           = nc_varget(undamped.fname,'VS');
undamped.VS    = [myVS(:,73:144)  myVS(:,1:73)]; % swap the hemispheres, repeat Greenwich


% damped ... 

damped.fname = '/gpfs/ptmp/raeder/Single_CAM/Halo-diff_12-5-2000/caminput_VS_266.nc';
damped.lon   = nc_varget(damped.fname,'slon');
damped.lat   = nc_varget(damped.fname,'lat');
damped.lev   = nc_varget(damped.fname,'lev');
damped.time  = nc_varget(damped.fname,'time');
myVS         = nc_varget(damped.fname,'VS');
damped.VS    = [myVS(:,73:144)  myVS(:,1:73)]; % swap the hemispheres

% difference ... 

diffed.fname = '/gpfs/ptmp/raeder/Single_CAM/Halo-diff_12-5-2000/1098_Prior-caminput.nc';
diffed.lon   = nc_varget(diffed.fname,'slon');
diffed.lat   = nc_varget(diffed.fname,'lat');
diffed.lev   = nc_varget(diffed.fname,'lev');
diffed.time  = nc_varget(diffed.fname,'time');
myVS         = nc_varget(diffed.fname,'VS');
diffed.VS    = [myVS(:,73:144)  myVS(:,1:73)]; % swap the hemispheres

%% Make some plots
% To be comparable with the other figures in the article, the hemispheres 
% need to be [-180,180], the major ticks and labels are every 30 degrees.
% no titles.
% bigger fonts
% Peter wanted .eps or .ps - good - can use 'painters'

datmat = [undamped.VS(:); damped.VS(:)];
datmin = min(datmat);
datmax = max(datmat);
clim = [datmin datmax];

% to reflect the swap
dlon = mean(diff(undamped.lon));
slon = [undamped.lon(73)-180:dlon: undamped.lon(144)+2*dlon] - 180;
slon = [-180 180];

% define matching labels

xticks = [-180:30:180];
yticks = [ -90:30:90 ];

xticklabels = ['180 ';'150W';'120W';'90W ';'60W ';'30W ';'  0 ';...
               '30E ';'60E ';'90E ';'120E';'150E';'180 '];
yticklabels = ['90S';'60S';'30S';' 0 ';'30N';'60N';'90N'];

figure(1); clf; orient landscape; colormap(bob); set(gcf,'renderer','painters')
imagesc(slon, undamped.lat, undamped.VS, clim)
set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
set(gca,'XTick',xticks,'XTickLabel',xticklabels)
set(gca,'YTick',yticks,'YTickLabel',yticklabels)
% title('Meridional wind with gridpoint-scale noise')
worldmap('hollow','dateline');
axis image
h = colorbar('vert');
set(h,'FontSize',14)
set(get(h,'YLabel'),'String','m/s','FontSize',14);

figure(2); clf; orient landscape; colormap(bob); set(gcf,'renderer','painters')
imagesc(slon, damped.lat, damped.VS, clim)
set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
set(gca,'XTick',xticks,'XTickLabel',xticklabels)
set(gca,'YTick',yticks,'YTickLabel',yticklabels)
% title('Meridional wind with del^4 damping')
worldmap('hollow','dateline');
axis image
h = colorbar('vert');
set(h,'FontSize',14)
set(get(h,'YLabel'),'String','m/s','FontSize',14);

figure(3); clf; orient landscape; colormap(bob); set(gcf,'renderer','painters')
imagesc(slon, diffed.lat, diffed.VS)
set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
set(gca,'XTick',xticks,'XTickLabel',xticklabels)
set(gca,'YTick',yticks,'YTickLabel',yticklabels)
% title('difference')
worldmap('hollow','dateline');
axis image
h = colorbar('vert');
set(h,'FontSize',14)
set(get(h,'YLabel'),'String','m/s','FontSize',14);

% print the graphics

print(1,'-depsc','VS_undamped')
print(2,'-depsc','VS_damped')
print(3,'-depsc','undamped-damped')

% For whatever reason, the 'painters' renderer is just miserable with
% pdf files - but great with postscript. Both OpenGL and zbuffer seem
% to generate decent pdf files - the default font size is a bit fuzzy,
% so I'm not sure these aren't being rasterized.

set(1,'renderer','OpenGL'); print(1,'-dpdf','VS_undamped');
set(2,'renderer','OpenGL'); print(2,'-dpdf','VS_damped');
set(3,'renderer','OpenGL'); print(3,'-dpdf','undamped-damped');

% Turns out - index 113 (of 256) is precisely zero on the difference plot.
% figure(3)
% bob = gauss3; bob(113,:) = 1; colormap(bob)
% bob = gauss3; bob(112:114,:) = 1; colormap(bob)
% bob = gauss3; bob(111:115,:) = 1; colormap(bob)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
