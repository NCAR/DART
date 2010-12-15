% PeterFigures generates 'publication-quality' figures for Peter Lauritzen
%
%

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

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
% is not what we have. Should swap halves to be consistent

% Undamped ... the noisy one.

undamped.fname = '/gpfs/ptmp/raeder/Single_CAM/Halo-diff_12-5-2000/Prior_memb1_VS_266.nc';
undamped.lon   = nc_varget(undamped.fname,'slon');
undamped.lat   = nc_varget(undamped.fname,'lat');
undamped.lev   = nc_varget(undamped.fname,'lev');
undamped.time  = nc_varget(undamped.fname,'time');
undamped.VS    = nc_varget(undamped.fname,'VS');

% damped ... 

damped.fname = '/gpfs/ptmp/raeder/Single_CAM/Halo-diff_12-5-2000/caminput_VS_266.nc';
damped.lon   = nc_varget(damped.fname,'slon');
damped.lat   = nc_varget(damped.fname,'lat');
damped.lev   = nc_varget(damped.fname,'lev');
damped.time  = nc_varget(damped.fname,'time');
damped.VS    = nc_varget(damped.fname,'VS');

% difference ... 

diffed.fname = '/gpfs/ptmp/raeder/Single_CAM/Halo-diff_12-5-2000/1098_Prior-caminput.nc';
diffed.lon   = nc_varget(diffed.fname,'slon');
diffed.lat   = nc_varget(diffed.fname,'lat');
diffed.lev   = nc_varget(diffed.fname,'lev');
diffed.time  = nc_varget(diffed.fname,'time');
diffed.VS    = nc_varget(diffed.fname,'VS');

%% Make some plots

datmat = [undamped.VS(:); damped.VS(:)];
datmin = min(datmat);
datmax = max(datmat);
clim = [datmin datmax];

figure(1); clf; orient landscape; colormap(bob)
imagesc(undamped.lon, undamped.lat, undamped.VS, clim)
set(gca,'YDir','normal','TickDir','out','XMinorTick','on')
title('Meridional wind with gridpoint-scale noise')
worldmap;
colorbar('vert')

figure(2); clf; orient landscape; colormap(bob)
imagesc(damped.lon, damped.lat, damped.VS, clim)
set(gca,'YDir','normal','TickDir','out','XMinorTick','on')
title('Meridional wind with del^4 damping')
worldmap;
colorbar('vert')

figure(3); clf; orient landscape; colormap(bob)
imagesc(diffed.lon, diffed.lat, diffed.VS)
set(gca,'YDir','normal','TickDir','out','XMinorTick','on')
title('difference')
worldmap;
colorbar('vert')

print(1,'-dpng','VS_undamped')
print(2,'-dpng','VS_damped')
print(3,'-dpng','undamped-damped')

print(1,'-dpsc','VS_undamped')
print(2,'-dpsc','VS_damped')
print(3,'-dpsc','undamped-damped')

% Turns out - index 113 (of 256) is precisely zero on the difference plot.
% figure(3)
% bob = gauss3; bob(113,:) = 1; colormap(bob)
% bob = gauss3; bob(112:114,:) = 1; colormap(bob)
% bob = gauss3; bob(111:115,:) = 1; colormap(bob)
