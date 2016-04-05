function obs = plot_coverage(fname)
%% plot_coverage examines the spatial and temporal coverage for an observation type.
%  obs_seq_coverage must create the netCDF file read by this routine.
% 
% fname = 'obs_coverage.nc';
%
% plot_coverage(fname)

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

locs    = nc_varget(fname,'location');

obs.fname         = fname;
obs.ObsTypeString = nc_read_att(fname,nc_global,'obs_of_interest');
obs.stations      = nc_varget(fname,'stations');
obs.ntimes        = nc_varget(fname,'ntimes');
obs.times         = nc_varget(fname,'time');
obs.timeunits     = nc_read_att(fname,'time','units');

obs.lons       = locs(:,1);
obs.lats       = locs(:,2);
obs.vert       = locs(:,3);
obs.region     = [ min(obs.lons) max(obs.lons) min(obs.lats) max(obs.lats) ];
timebase       = sscanf(obs.timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
obs.timeorigin = datenum(timebase(1),timebase(2),timebase(3));
obs.times      = obs.times + obs.timeorigin;

allt     = obs.times(:);
obs.tmin = datestr(min(allt));
obs.tmax = datestr(max(allt));
clear allt;

inds = find(obs.ntimes == 8);
fprintf('There are %d stations with precisely 8 observation times.\n',length(inds))

%% Now pack the data in the same fashion as the cell array of column labels.

obs.lonindex  = 1;
obs.latindex  = 2;
obs.zindex    = 3;
obs.ntindex   = 4;
obs.Iindex    = 5;

obs.colnames{obs.lonindex}  = 'longitude';
obs.colnames{obs.latindex}  = 'latitude';
obs.colnames{obs.zindex}    = 'vertical';
obs.colnames{obs.ntindex}   = 'num times';
obs.colnames{obs.Iindex}    = 'stationID';

global obsmat
obsmat = zeros(length(obs.lons),5);
obsmat(:,obs.lonindex ) = obs.lons;
obsmat(:,obs.latindex ) = obs.lats;
obsmat(:,obs.zindex   ) = obs.vert;
obsmat(:,obs.ntindex  ) = obs.ntimes;
obsmat(:,obs.Iindex   ) = obs.stations;

global timemat
timemat = obs.times;

%% create the linked plots

figure1 = figure(1); clf(figure1); orient tall; wysiwyg

%----------------------------------------------------------------------
% Top half is a 2D scatterplot
%----------------------------------------------------------------------

axes0 = axes('Parent',figure1,'OuterPosition',[0 0.4 1 0.57],'FontSize',16);
view(axes0,[0 90]);
grid(axes0,'on');
hold(axes0,'all');

xstring = sprintf('obsmat(:,%d)',obs.lonindex);
ystring = sprintf('obsmat(:,%d)',obs.latindex);

scatter(obsmat(:,obs.lonindex), obsmat(:,obs.latindex), 18, obsmat(:,obs.ntindex), 'filled', ...
             'Parent',axes0,'DisplayName','observation locations', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
% myworldmap(obs);

xlabel(obs.colnames{obs.lonindex});
ylabel(obs.colnames{obs.latindex});

hc = colorbar;
set(get(hc,'YLabel'),'String','# of obs times','Interpreter','none')

ntmax = length(obs.times);
h = title({sprintf('%s locations',obs.ObsTypeString), ...
           sprintf('%s ---> %s',obs.tmin,obs.tmax) });
set(h,'Interpreter','none')

%----------------------------------------------------------------------
% Bottom half is a "histogram"
%----------------------------------------------------------------------

axes1 = axes('Parent',figure1,'OuterPosition',[0 0.0 1 0.38],'FontSize',14);
box(axes1,'on');
hold(axes1,'all');

xstring = sprintf('obsmat(:,%d)',obs.Iindex);
ystring = sprintf('obsmat(:,%d)',obs.ntindex);

scatter(obsmat(:,obs.Iindex), obsmat(:,obs.ntindex), ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
xlabel('''station''')
ylabel('# of obs times')

refreshdata
linkdata on

%% Create figure for ancillary plots
%  ====================================================================

figure2 = figure(2); clf(figure2); orient tall; wysiwyg

%% Create axes for time VS. station
axes2 = axes('Parent',figure2,'OuterPosition',[0 0 1 1],'FontSize',14);
set(axes2,'XAxisLocation','bottom')
box(axes2,'on');
hold(axes2,'all');

hist(obsmat(:,obs.ntindex),single([1:max(obsmat(:,obs.ntindex))]))
h = title(obs.ObsTypeString);
set(h,'Interpreter','none')
xlabel('number of observations at a location')
ylabel('count')

%% Create figure for ancillary plots
%  ====================================================================

figure3 = figure(3); clf(figure3); orient tall; wysiwyg

%% Create axes for time VS. station
axes3 = axes('Parent',figure3,'OuterPosition',[0 0 1 1],'FontSize',14);
set(axes3,'XAxisLocation','bottom')
box(axes3,'on');
hold(axes3,'all');

ystring = sprintf('obsmat(:,%d)',obs.Iindex);

h = plot(timemat, obsmat(:,obs.Iindex),'kd','Parent',axes3, ...
             'DisplayName','time vs station', ...
             'YDataSource',ystring, ...
             'XDataSource','timemat');

ax = axis;
datetick(axes3,'x',20);
axis([min(timemat(:))-0.25 max(timemat(:))+0.25 ax(3) ax(4)])
h = title(obs.ObsTypeString);
set(h,'Interpreter','none')
xlabel('Date - DD/MM/YY')
ylabel('station ID')

%% thats it folks

% refreshdata
% linkdata on



function h = myworldmap(obs)

%%--------------------------------------------------------------------------
% GET THE ELEVATION DATA AND SET UP THE ASSOCIATED COORDINATE DATA
%---------------------------------------------------------------------------

load topo;               % GET Matlab-native [180x360] ELEVATION DATASET
lats = -89.5:89.5;       % CREATE LAT ARRAY FOR TOPO MATRIX
lons = 0.5:359.5;        % CREATE LON ARRAY FOR TOPO MATRIX
nlon = length(lons);
nlat = length(lats);

%%--------------------------------------------------------------------------
% IF WE NEED TO SWAP HEMISPHERES, DO SO NOW.
% If we didn't explicitly tell it, make a guess.
%---------------------------------------------------------------------------

ax   = axis;

if (ax(1) < -2)
   lons = lons - 180.0;
   topo = [ topo(:,nlon/2+1:nlon) topo(:,1:nlon/2) ];
end

%%--------------------------------------------------------------------------
% We need to determine the geographic subset of the elevation matrix.
%---------------------------------------------------------------------------

lon_ind1 = find(ax(1) <= lons, 1);
lon_ind2 = find(ax(2) <= lons, 1);
lat_ind1 = find(ax(3) <= lats, 1);
lat_ind2 = find(ax(4) <= lats, 1);

if (isempty(lon_ind1)), lon_ind1 = 1;    end;
if (isempty(lon_ind2)), lon_ind2 = nlon; end;
if (isempty(lat_ind1)), lat_ind1 = 1;    end;
if (isempty(lat_ind2)), lat_ind2 = nlat; end;

elev = topo(lat_ind1:lat_ind2,lon_ind1:lon_ind2);
x    = lons(lon_ind1:lon_ind2);
y    = lats(lat_ind1:lat_ind2);

%%--------------------------------------------------------------------------
% Contour the "subset" - and give the whole thing an appropriate zlevel
% so the continents are either at the top of the plot (for depth), or
% the bottom (for just about everything else.
% There are differences between 6.5 and 7.0 that make changing the colors
% of the filled contours a real pain.
%---------------------------------------------------------------------------

orgholdstate = ishold;
hold on;

switch  lower(obs.Zunits)
   case 'depth'
      set(gca,'Zdir','reverse')
      zlevel = min(ax(5:6));
   case 'pressure'
      zlevel = max(ax(5:6));
      set(gca,'Zdir','reverse')
   otherwise
      set(gca,'Zdir','normal')
      zlevel = min(ax(5:6));
end

fcolor = [0.7 0.7 0.7];    % light grey

[c,h] = contourf(x,y,elev,[0.0 0.0],'k-');

h_patch = get(h, 'Children');

for i = 1:numel(h_patch)
    y = get(h_patch(i), 'YData');
    s = size(y);
    set(h_patch(i), 'ZData', zlevel*ones(s),'FaceColor',fcolor);
    set(h_patch(i),'AlphaDataMapping','none','FaceVertexAlphaData',0.3)
    set(h_patch(i),'FaceAlpha',0.3)
end

if (orgholdstate == 0), hold off; end;


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

