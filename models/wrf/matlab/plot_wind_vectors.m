function obs = plot_wind_vectors(fname,ncname,platform,levelindex)
% plot_wind_vectors(fname,ncname,platform,ilevel)
%
% fname    = 'wind_vectors.006.dat';
% ncname   = 'obs_diag_output.nc';
% platform = 'SAT';
% ilevel   = 1;         % surface == 1?
% ilevel   = -1;        % all levels == -1
%
% obs = plot_wind_vectors(fname,ncname,platform,ilevel);

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

scalefactor = 10.0;

f = netcdf(ncname,'nowrite');
platforms = f{'ObservationTypes'}(:);
timeunits = f{'time'}.units(:);
close(f);

timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
obs        = Load_Parse(fname, platforms, platform, levelindex );

if  isempty(obs)
   clf;
   axis([0 1 0 1]); axis image
   h = text(0.5,0.5,sprintf('%s has no %s data',fname,platform));
   set(h,'Interpreter','none')
   return
end

obs.times  = obs.times + timeorigin;
t1         = datestr(min(obs.times));
t2         = datestr(max(obs.times));

clf;
axlims     = DrawBackground(obs,scalefactor);

goodUV = find( (obs.Uqc < 1) & (obs.Vqc < 1));
baadUV = find( (obs.Uqc > 1) & (obs.Vqc > 1));
goodU  = find( (obs.Uqc < 1) & (obs.Vqc > 0));
goodV  = find( (obs.Vqc < 1) & (obs.Uqc > 0));

legh   = [];
legstr = {};

if ~ isempty(goodUV)
   hgood = obs2plot(obs,goodUV,[0 0 0],scalefactor);
   legh   = [legh; hgood];
   legstr{length(legstr)+1} = sprintf('%d good',length(goodUV));
end

if ~ isempty(baadUV)
   hbaadUV = obs2plot(obs,baadUV,[1 0 0],scalefactor);
   legh   = [legh; hbaadUV];
   legstr{length(legstr)+1} = sprintf('%d badUbadV',length(baadUV));
end

if ~ isempty(goodU)
   hgoodU = obs2plot(obs,goodU,[0 1 0],scalefactor);
   legh   = [legh; hgoodU];
   legstr{length(legstr)+1} = sprintf('%d goodUbadV',length(goodU));
end

if ~ isempty(goodV)
   hgoodV = obs2plot(obs,goodV,[0 0 1],scalefactor);
   legh   = [legh; hgoodV];
   legstr{length(legstr)+1} = sprintf('%d badUgoodV',length(goodV));
end

h = title({sprintf('%s %s %s',t1,platform,t2), ...
           sprintf('binlevel %d',levelindex)});
set(h,'FontSize',18)

h = xlabel(fname); set(h,'Interpreter','none');

legend(legh,legstr,'Location','NorthWestOutside')

hold off;

function axlims = DrawBackground(obs, xfactor)
%----------------------------------------------------------------------

% Figure out bounds of the data
axlims = [ min(obs.lon) max(obs.lon) min(obs.lat) max(obs.lat) ];

% It is nice to have a little padding around the perimeter
dx = 0.05 * (axlims(2) - axlims(1));
dy = 0.05 * (axlims(4) - axlims(3));

% axis(axlims + [-dx dx -dy dy])
axis([115 145 5 35])
axis image
axlims = axis;

% It is nice to know where the land is
worldmap('light');
hold on;

tx = axlims(1)+2;
ty = axlims(4)-2;

U = 10/xfactor; V =  0/xfactor;
h = quiver(tx, ty, U, V, 0.0 ,'LineWidth',4.0,'Color','k');
set(h,'LineWidth',3.0)

U =  0/xfactor; V = 10/xfactor;
h = quiver(tx, ty, U, V, 0.0 ,'LineWidth',4.0,'Color','k');
set(h,'LineWidth',3.0)
h = text(tx, ty,'10 m/s','VerticalAlignment','top');



function h1 = obs2plot(obs, mask, colspec, xfactor)
%----------------------------------------------------------------------
lon    = obs.lon(mask);
lat    = obs.lat(mask);
U      = obs.U(mask)/xfactor;
V      = obs.V(mask)/xfactor;
h1     = quiver(lon, lat, U, V, 0.0);
h2     = plot(lon, lat, '.','MarkerSize',4);
set(h1,'Color',colspec)
set(h2,'Color',colspec)



function newobs = Load_Parse(fname, platforms, platform, levelindex)
%----------------------------------------------------------------------
% Makes no attempt to find/replace/identify MISSING values

obsmat = load(fname);

platformIDs = unique(obsmat(:,1));
uid         = floor(platformIDs/100);
vid         = platformIDs - uid*100;
Ustrings    = platforms(uid,:);
Vstrings    = platforms(vid,:);

nplatforms = length(uid);
pid        = [];
newobs     = [];

for i = 1:nplatforms
   uindex = findstr(Ustrings(i,:),'_U_WIND_COMPONENT') - 1;
   vindex = findstr(Vstrings(i,:),'_V_WIND_COMPONENT') - 1;

   Ubase  = Ustrings(i,1:uindex);
   Vbase  = Vstrings(i,1:vindex);

   if ~strcmp(Ubase,Vbase)
      error('U and V wind component platforms do not match')
   end

   if strcmp(lower(Ubase),lower(platform))
      pid = platformIDs(i);
   end

   inds = find(obsmat(:,1) == platformIDs(i));
   nobs = length(inds);

   disp(sprintf('%6d %14s observations in %s (%4d)',nobs,Ubase,fname,platformIDs(i)))

end

if isempty(pid)
   disp(sprintf('no %s observations in %s',platform,fname))
   return
end

inds   = find(obsmat(:,1) == pid);

if isempty(inds)
   disp(sprintf('no observations of type %d in %s',pid,fname))
end

obs.platform = obsmat(inds, 1);
obs.day      = obsmat(inds, 2);
obs.seconds  = obsmat(inds, 3);
obs.lon      = obsmat(inds, 4);
obs.lat      = obsmat(inds, 5);
obs.level    = obsmat(inds, 6);
obs.levind   = obsmat(inds, 7);
obs.Uqc      = obsmat(inds, 8);
obs.Vqc      = obsmat(inds, 9);
obs.U        = obsmat(inds,10);
obs.V        = obsmat(inds,11);
obs.Upr      = obsmat(inds,12);
obs.Vpr      = obsmat(inds,13);
obs.Upo      = obsmat(inds,14);
obs.Vpo      = obsmat(inds,15);
obs.times = obs.day + obs.seconds/86400;

% Sort by level

if ( levelindex > 0 ) 
   levelinds = find(obs.levind == levelindex);
else
   levelinds = 1:length(inds);
end


if isempty(levelinds)
   error('no observations of type %d at level %d',pid,levelindex)
else
   disp(sprintf('%d observations of type %d at level %d',length(levelinds),pid,levelindex))
end

clear newobs

newobs.platform = obs.platform(levelinds);
newobs.day      = obs.day(     levelinds);
newobs.seconds  = obs.seconds( levelinds);
newobs.lon      = obs.lon(     levelinds);
newobs.lat      = obs.lat(     levelinds);
newobs.level    = obs.level(   levelinds);
newobs.levind   = obs.levind(  levelinds);
newobs.Uqc      = obs.Uqc(     levelinds);
newobs.Vqc      = obs.Vqc(     levelinds);
newobs.U        = obs.U(       levelinds);
newobs.V        = obs.V(       levelinds);
newobs.Upr      = obs.Upr(     levelinds);
newobs.Vpr      = obs.Vpr(     levelinds);
newobs.Upo      = obs.Upo(     levelinds);
newobs.Vpo      = obs.Vpo(     levelinds);
newobs.times    = obs.times(   levelinds);

clear obsmat obs
