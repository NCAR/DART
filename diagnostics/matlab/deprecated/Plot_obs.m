%% Plot_obs
%
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

rad2deg = 45/atan(1);

QTY_U   =  1;
QTY_V   =  2;
QTY_PS  =  3;
QTY_T   =  4;
QTY_VR  = 12;
QTY_REF = 13;
QTY_U10 = 14;
QTY_V10 = 15;
QTY_T2  = 16;
QTY_Q2  = 17;
QTY_TD2 = 18;

var = input('Input obs to plot (VR, REF): ');
if strcmp(var,'REF')
   var_units = 'dBZ';
end
if strcmp(var,'VR')
   var_units = 'm/s';
end
lev = input('Input lev (m): ');
tol = input('Input margin (m): ');

map_proj = {'lambert', 'ups', 'mercator'};

if ( exist('preassim.nc','file') ~= 0 )
  fname = 'preassim';
elseif ( exist('postassim.nc','file') ~= 0 )
  fname = 'postassim';
elseif ( exist('true_state.nc','file') ~= 0 )
  fname = 'true_state';
else
  fname = input('Enter the name of the netCDF file containing domain information: ');
end

if strcmp(fname(1:8),'wrfinput')

     nc = netcdf( fname , 'nowrite' ) ;

     stdlat1 = nc.TRUELAT1(:);
     stdlat2 = nc.TRUELAT2(:);
     cen_lat = nc.CEN_LAT(:);
     cen_lon = nc.CEN_LON(:);
     mp = nc.MAP_PROJ(:);

     close(nc);

     id = 1;

     xlon = ncread(fname, 'XLONG');
     xlat = ncread(fname, 'XLAT');

else

     stdlat1 = ncread(fname, 'TRUELAT1');
     stdlat2 = ncread(fname, 'TRUELAT2');
     cen_lat = ncread(fname, 'CEN_LAT');
     cen_lon = ncread(fname, 'CEN_LON');
     mp      = ncread(fname, 'MAP_PROJ');

     num_domains = size(mp,1);

if (num_domains > 1)

   disp(['Number of domains: ',int2str(num_domains)])
   id = input('Input domain id: ');

else

   id = 1;

end

     xlon = ncread(fname, ['XLON_d0',int2str(id)]);
     xlat = ncread(fname, ['XLAT_d0',int2str(id)]);

end

we = size(xlon, 2);
sn = size(xlat, 1);
minlat = min(xlat(:)); maxlat = max(xlat(:));
minlon = min(xlon(:)); maxlon = max(xlon(:));

if(exist('a','var') == 0)
   a = ReadASCIIObsSeq('obs_seq.final');
end

days = -1;
secs = -1;
j = 0;
disp('Available times')
disp('day     second')
for i = 1:a.num_obs
   if (a.days(i) ~= days | a.secs(i) ~= secs)
     j = j + 1;
     day(j) = a.days(i);
     sec(j) = a.secs(i);
     beg(j) = i;
     days = a.days(i);
     secs = a.secs(i);
     disp([int2str(days),' ',int2str(secs),' ',int2str(j)])
   end
end
beg(j+1) = a.num_obs + 1;

j = input(['Enter period # of interest ( 1 - ',int2str(j),' ): ']);

m = min([2 a.num_copies]);

string1 = sprintf('%s (%s)  %d (m)  day = %d  sec = %d  ',var,var_units,lev,a.days(beg(j)),a.secs(beg(j)));

metadata = {'observations', 'prior ensemble mean'};

for pane = 1:m,

iu=0;
iv=0;
it=0;
iu10=0;
iv10=0;
it2=0;
itd2=0;
ips=0;
ivr=0;
iref=0;

field = zeros(sn,we);

field(:,:) = NaN ;

figure(pane);
%subplot(1,m,pane);

axesm(map_proj{mp(id)},'Origin',[0 cen_lon(id) 0],'MapParallels',[stdlat1(id) stdlat2(id)],...
      'MapLatLimit',[minlat maxlat],'MapLonLimit',[minlon maxlon]); framem;

[xlim ylim]=mfwdtran([xlat(1,1) xlat(sn,we)],[xlon(1,1) xlon(sn,we)]);
set(gca,'xlim',[min(xlim(:)) max(xlim(:))]);
set(gca,'ylim',[min(ylim(:)) max(ylim(:))]);

plotm(coast,'color',[0 0 0]);
plotm(usalo('statebvec'),'color',[0 0 0]);
plotm(usalo('conusvec'),'color',[0 0 0]);

%title({string1,a.metadata(pane)}, 'Fontsize',12)
title([string1,metadata{pane}], 'Fontsize',12)

for i = beg(j):min([(beg(j+1)-1) a.num_obs]),
%for i = 1:180,

   lon = rad2deg*a.loc(i,1);
   lat = rad2deg*a.loc(i,2);
   hei = a.loc(i,3);

   if ( a.kind(i) == QTY_U )

     iu = iu + 1;

   scatterm(lat,lon,'xb')

   elseif ( a.kind(i) == QTY_V )

     iv = iv + 1;

   elseif ( a.kind(i) == QTY_T )

     it = it + 1;

   scatterm(lat,lon,'+r')

   elseif ( a.kind(i) == QTY_U10 )

     iu10 = iu10 + 1;

     uwind = a.obs(pane,i);

   elseif ( a.kind(i) == QTY_V10 )

     iv10 = iv10 + 1;

     if ( (uwind ~= -888888.0) && (a.obs(pane,i) ~= -888888.0) )
%        quiverm(lat,lon,a.obs(pane,i),uwind)
     end

   elseif ( a.kind(i) == QTY_T2 )

     if (a.obs(pane,i) ~= -888888.0)
       textm(lat,lon,[num2str(round(a.obs(pane,i)-273.15)),' '],...
     'VerticalAlignment','bottom',...
     'HorizontalAlignment','right')
     end

     it2 = it2 + 1;

   scatterm(lat,lon,'og')

   elseif ( a.kind(i) == QTY_TD2 )

     if (a.obs(pane,i) ~= -888888.0)
       textm(lat,lon,[num2str(round(a.obs(pane,i)-273.15)),' '],...
     'VerticalAlignment','top',...
     'HorizontalAlignment','right')
     end

     itd2 = itd2 + 1;

   scatterm(lat,lon,'og')

   elseif ( a.kind(i) == QTY_PS )

     if (a.obs(pane,i) ~= -888888.0)
       textm(lat,lon,[' ',num2str(round(a.obs(pane,i)/100))],...
     'VerticalAlignment','bottom',...
     'HorizontalAlignment','left')
     end

     ips = ips + 1;

   scatterm(lat,lon,'og')

   elseif ( a.kind(i) == QTY_VR )

     ivr = ivr + 1;

     if strcmp(var,'VR')
        if((hei > (lev-tol)) && (hei < (lev+tol)))
          scatterm(lat,lon,16,a.obs(pane,i),'s','filled')
        end
     end

   elseif ( a.kind(i) == QTY_REF )

     iref = iref + 1;

     if strcmp(var,'REF')
        if(hei > (lev-tol) && hei < (lev+tol))
     if(lon > 180)
        lon = lon - 360;
     end
     [C k] = min(reshape((xlat - lat).^2 + (xlon - lon).^2,we*sn,1));
%     field(k) = a.obs(pane,i);
     if a.obs(pane,i) ~= 0.0
        field(k) = 10.0*log10(a.obs(pane,i));
     else
     field(k) = -16.0;
     end
        end
     end

   else

     disp(['Kind for obs ',int2str(i),' is ',int2str(a.kind(i))])

   end

   hold on

end

if min(min(field)) ~= max(max(field))
   if strcmp(var,'REF')
     eval('dbz_colors')
   end
   h = pcolorm(xlat,xlon,field);
%   caxis([min(iso(:)),max(iso(:))]);
end

cb = colorbar('vert'); set(cb,'Fontsize',12);

   wysiwyg

fprintf('# of U   %d\n', iu)
fprintf('# of V   %d\n', iv)
fprintf('# of T   %d\n', it)
fprintf('# of Vr  %d\n', ivr)
fprintf('# of Ref %d\n', iref)
fprintf('# of U10 %d\n', iu10)
fprintf('# of V10 %d\n', iv10)
fprintf('# of T2  %d\n', it2)
fprintf('# of TD2 %d\n', itd2)
fprintf('# of PS  %d\n', ips)

end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
