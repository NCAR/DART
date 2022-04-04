function CAM_DART_correl(datadir,DARTfile,DARTvarname,DARTlevel,CAMvarname,CAMlocation)
%% CAM_DART_correl  explores correlation between one CAM location and the DART ensemble
%
% datadir     = '/ptmp/raeder/Cam3.6/Hist0';
% DARTfile    = 'Pr_06_ens.nc';
% DARTvarname = 'T';
% DARTlevel   = 600;
% CAMvarname  = 'EVAPPREC';
% CAMlocation = [170, 31, 600];  % lon/lat/level
%
% CAM_DART_correl(datadir,DARTfile,DARTvarname,DARTlevel,CAMvarname,CAMlocation)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

posmat = [0.1 0.60 0.8 0.25;
          0.1 0.10 0.8 0.40];

%----------------------------------------------------------------------
% Get the time series from the CAM files:
% EVAPPREC(time, lev, lat, lon) ;
%----------------------------------------------------------------------

% For this variable, should query for coordinate variables
% not handling variables on staggered grids correctly as it stands.

camname = sprintf('%s/FV_2deg-noleap-O2-Dev20-1-%02d.cam2.h0.2003-07-08-21600_%02d.nc',datadir,1,1);
CAMlevels = nc_varget(camname,'lev'); CAMNlevels = length(CAMlevels);
CAMlats   = nc_varget(camname,'lat'); CAMNlats   = length(CAMlats);
CAMlons   = nc_varget(camname,'lon'); CAMNlons   = length(CAMlons);
CAMunits  = nc_attget(camname,CAMvarname,'units'); 
CAMtimestr = timebase(camname);

diffs  = abs(CAMlevels - CAMlocation(3)); [camlevel,camlevind] = min(diffs);
diffs  = abs(CAMlats   - CAMlocation(2)); [camlat,  camlatind] = min(diffs);
diffs  = abs(CAMlons   - CAMlocation(1)); [camlon,  camlonind] = min(diffs);

camlevel = CAMlevels(camlevind);
camlat   = CAMlats(  camlatind);
camlon   = CAMlons(  camlonind);

% start = [ 0 camlevind-1 camlatind-1 camlonind-1];
% count = [-1     1           1           1      ];

start = [ 0 camlevind-1   0   0];
count = [ 1     1        -1  -1];

x = zeros(80,1);
cammean = zeros(CAMNlats,CAMNlons);

for ensmem = 1:80

   camname = sprintf('%s/FV_2deg-noleap-O2-Dev20-1-%02d.cam2.h0.2003-07-08-21600_%02d.nc',datadir,ensmem,ensmem);

   bob       = nc_varget(camname,CAMvarname,start,count);
   cammean   = cammean + bob;
   x(ensmem) = bob(camlatind,camlonind); 

end

cammean = cammean/80;

% should check to make sure they are no 'MISSING' values in the x array ...

%----------------------------------------------------------------------
% Get the DART fields.
% float Q(time, copy, lat, lon, lev) ;
%----------------------------------------------------------------------

prfname    = sprintf('%s/%s',datadir,DARTfile);
DARTlons   = nc_varget(prfname,'lon');
DARTlats   = nc_varget(prfname,'lat');
DARTlevels = nc_varget(prfname,'lev');

diffs = abs(DARTlevels - DARTlevel);
[dartlevel,dartlevind] = min(diffs);
dartlevel = DARTlevels(dartlevind);

copyind1 = get_copy_index(prfname,'ensemble member 1');
copyindN = get_copy_index(prfname,'ensemble member 80');
copyN = copyindN - copyind1 + 1;

start = [ 0  copyind1-1  0  0 dartlevind-1];
count = [-1  copyN      -1 -1      1];
y     = nc_varget(prfname,DARTvarname,start,count);

[Nmem, Nlat, Nlon] = size(y);

%----------------------------------------------------------------------
% Reshape the DART data to be in state-space form ... Nt-by-Nx
% perform the correlation,
% reshape the correlation coefficients back to a Lat-Lon matrix.
%----------------------------------------------------------------------

Yss = reshape(y,[Nmem Nlat*Nlon]); 
a   = corr(x,Yss);
b   = reshape(a,[Nlat Nlon]);

clear Yss a

bob = b(:);
abob = abs(bob);
datamax = max(abob);

%datamax = 0.5;
%----------------------------------------------------------------------
% make the plots
%----------------------------------------------------------------------

figure(1); clf; orient tall

str1 = sprintf('%s correlated against %s',CAMvarname,DARTvarname);
str2 = sprintf('%s [%.2f , %.2f] level %.2f',CAMvarname,camlon,camlat,camlevel);
str3 = sprintf('%s level %.2f',DARTvarname,dartlevel);

subplot('position',posmat(1,:))
plot(x)
title( {str2, CAMtimestr} )
ylabel(CAMunits)
xlabel('Ensemble Member');

subplot('position',posmat(2,:))
mymap = [1 1 1; jet];
colormap(mymap)
imagesc(CAMlons,CAMlats,cammean)
set(gca,'YDir','normal')
h = colorbar;
set(get(h,'YLabel'),'String',CAMunits)
worldmap;
axis image
hold on; 
plot([  0 360],[camlat camlat],'k', ...
     [camlon camlon],[-90 90],'k');
title(sprintf('Ensemble mean %s level %.2f %s',CAMvarname,camlevel,CAMtimestr))


figure(2); clf; orient tall


subplot('position',posmat(2,:))
imagesc(DARTlons,DARTlats,b,[-datamax datamax])
set(gca,'YDir','normal')
mymap = jet;
mymap(30:35,:) = 1;
colormap(mymap);
colorbar
worldmap;
axis image
hold on; 
plot(camlon,camlat,'kd','MarkerSize',20,'LineWidth',2);

title({str1,str2,str3,CAMtimestr})




function timestring = timebase(fname)

t      = nc_varget(fname,'time');
units  = nc_attget(fname,'time','units'); 
mycal  = nc_attget(fname,'time','calendar'); 

timebase = sscanf(units,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
timestring = datestr(t + timeorigin);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
