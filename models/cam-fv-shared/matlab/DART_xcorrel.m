function DART_xcorrel(datadir,fname,varname1,location,varname2,level2)
%% DART_xcorrel  explores correlation between one DART location and the DART ensemble
%
% datadir  = '/ptmp/raeder/Cam3.6/Hist0';
% fname    = 'Pr_06_ens.nc';
% varname1 = 'T';
% location = [284, -54, 600];  % lon/lat/level
% varname2 = 'T';
% level2   = 600;
%
% DART_xcorrel(datadir,fname,varname1,location,varname2,level2)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

posmat = [0.1 0.60 0.8 0.25;
          0.1 0.10 0.8 0.40];

%----------------------------------------------------------------------
% Get the "time series" from the DART file:
%----------------------------------------------------------------------

% For this variable, should query for coordinate variables
% not handling variables on staggered grids correctly as it stands.

camname = sprintf('%s/%s',datadir,fname);

CAMlevels = nc_varget(camname,'lev'); CAMNlevels = length(CAMlevels);
CAMlats   = nc_varget(camname,'lat'); CAMNlats   = length(CAMlats);
CAMlons   = nc_varget(camname,'lon'); CAMNlons   = length(CAMlons);
CAMunits  = nc_attget(camname,varname1,'units'); 
CAMtimestr = timebase(camname);

diffs  = abs(CAMlevels - location(3)); [camlevel,camlevind] = min(diffs);
diffs  = abs(CAMlats   - location(2)); [camlat,  camlatind] = min(diffs);
diffs  = abs(CAMlons   - location(1)); [camlon,  camlonind] = min(diffs);

camlevel = CAMlevels(camlevind);
camlat   = CAMlats(  camlatind);
camlon   = CAMlons(  camlonind);

%   float T(time, copy, lat, lon, lev) ;
%                T:long_name = "Temperature" ;
%                T:units = "K" ;

copyind1 = get_copy_index(camname,'ensemble member 1');
copyindN = get_copy_index(camname,'ensemble member 80');
copyN = copyindN - copyind1 + 1;

start = [ 0  copyind1-1   0   0   camlevind-1];
count = [ 1  copyN       -1  -1        1     ];

bob = nc_varget(camname,varname1,start,count);
x   = bob(:,camlatind,camlonind); 

cammean = squeeze(mean(bob,1));

%----------------------------------------------------------------------
% Get the DART fields.
% float Q(time, copy, lat, lon, lev) ;
%----------------------------------------------------------------------

prfname    = sprintf('%s/%s',datadir,fname);
DARTlons   = nc_varget(prfname,'lon');
DARTlats   = nc_varget(prfname,'lat');
DARTlevels = nc_varget(prfname,'lev');

diffs = abs(DARTlevels - level2);
[dartlevel,dartlevind] = min(diffs);
dartlevel = DARTlevels(dartlevind);

start = [ 0  copyind1-1  0  0 dartlevind-1];
count = [ 1  copyN      -1 -1      1      ];
y     = nc_varget(prfname,varname2,start,count);

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

str1 = sprintf('%s correlated against %s',varname1,varname2);
str2 = sprintf('%s [%.2f , %.2f] level %.2f',varname1,camlon,camlat,camlevel);
str3 = sprintf('%s level %.2f',varname2,dartlevel);

subplot('position',posmat(1,:))
plot(x)
title( {str2, CAMtimestr} )
ylabel(CAMunits)
xlabel('Ensemble Member');

subplot('position',posmat(2,:))
mymap = jet;
colormap(mymap)
imagesc(CAMlons,CAMlats,cammean)
set(gca,'YDir','normal')
h = colorbar;
set(get(h,'YLabel'),'String',CAMunits)
worldmap;
axis image
hold on; 
plot([  0 360],[camlat camlat],'k', [camlon camlon],[-90 90],'k');
title(sprintf('Ensemble mean %s level %.2f   %s',varname1,camlevel,CAMtimestr))


figure(2); clf; orient tall


subplot('position',posmat(2,:))
imagesc(DARTlons,DARTlats,b,[-datamax datamax])
set(gca,'YDir','normal')
mymap = jet;
mymap(26:39,:) = 1;
colormap(mymap);
colorbar
worldmap;
axis image
hold on; 
plot([  0 360],[camlat camlat],'k', [camlon camlon],[-90 90],'k');

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
