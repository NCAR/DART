function linked_observations(obs)
%% linked_observations(obs) is a helper function for link_obs.m
% linked_observations is never meant to be called directly.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% obs has components (for example):
%            fname: '/ptmp/thoar/POP/CAM/POP8/obs_sequence_001.nc'
%    ObsTypeString: 'RADIOSONDE_TEMPERATURE'
%    ObsCopyString: 'NCEP BUFR observation'
%       CopyString: 'prior ensemble mean'
%         QCString: 'DART quality control'
%           region: [0 360 -90 90 -Inf Inf]
%          verbose: 1
%       timestring: [2x20 char]
%             lons: [794x1 double]
%             lats: [794x1 double]
%                z: [794x1 double]
%              obs: [794x1 double]
%             Ztyp: [794x1 int32]
%             keys: [794x1 int32]
%             time: [794x1 double]
%               qc: [794x1 int32]
%         colnames: {1x9 cell}
%         lonindex: 1
%         latindex: 2
%           zindex: 3
%         obsindex: 4
%        copyindex: 5
%          qcindex: 6
%         keyindex: 7
%        timeindex: 8
%         indindex: 9

global obsmat

%% ------------------------------------------------------------------------
%  Create figure and axes for 3D scatterplot
% -------------------------------------------------------------------------

figure1 = figure(1); clf(figure1);

fig1ax1 = axes('Parent',figure1,'OuterPosition',[0 0 1 0.90]);
view(fig1ax1,[-37.5 30]);
grid(fig1ax1,'on');

myworldmap(obs); hold on;

xstring = sprintf('obsmat(:,%d)',obs.lonindex);
ystring = sprintf('obsmat(:,%d)',obs.latindex);
zstring = sprintf('obsmat(:,%d)',obs.zindex  );

symbolsize = set_symbol_size(numel(obs.obs));

scatter3(obsmat(:,obs.lonindex), obsmat(:,obs.latindex), obsmat(:,obs.zindex), ...
             symbolsize, obsmat(:,obs.obsindex), 'filled', ...
             'Parent',fig1ax1, ...
             'DisplayName','observation locations', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring, ...
             'ZDataSource',zstring);

set(fig1ax1,'FontSize',18);

xlabel(obs.colnames{obs.lonindex});
ylabel(obs.colnames{obs.latindex});
zlabel(obs.colnames{obs.zindex});

tmin = min(obs.time);
tmax = max(obs.time);

h = title({obs.ObsTypeString, ...
      sprintf('"%s"',obs.ObsCopyString), ...
      sprintf('%s ---> %s',datestr(tmin),datestr(tmax)) });
set(h,'Interpreter','none')
linkdata on

%% ------------------------------------------------------------------------
% Create figure for ancillary plots
% -------------------------------------------------------------------------

figure2 = figure(2); clf(figure2); orient tall; wysiwyg

%% Create axes for time VS. QC
% FIXME ... choose a format for datetick based on the timespan
fig2ax1 = axes('Parent',figure2,'Position',[0.13 0.71 0.78 0.22]);
set(fig2ax1,'XAxisLocation','bottom')
box(fig2ax1,'on');

xstring = sprintf('obsmat(:,%d)',obs.timeindex);
ystring = sprintf('obsmat(:,%d)',obs.qcindex);
scatter(obsmat(:,obs.timeindex),obsmat(:,obs.qcindex), ...
    'Parent', fig2ax1, ...
    'DisplayName', 'time vs qc', ...
    'XDataSource', xstring, ...
    'YDataSource', ystring);

set(fig2ax1,'FontSize',14);

ax = axis;
ax(3:4) = [0 8]; %min/max range of DART QC values
axis(ax)
datetick(fig2ax1,'x')
ylabel(obs.colnames{obs.qcindex});
grid(fig2ax1,'on');

%% Create axes for observation count/density VS. time 

fig2ax2 = axes('Parent',figure2,'Position',[0.13 0.41 0.78 0.22]);
box(fig2ax2,'on');

xstring = sprintf('obsmat(:,%d)',obs.timeindex);
ystring = sprintf('obsmat(:,%d)',obs.indindex);
scatter(obsmat(:,obs.timeindex),obsmat(:,obs.indindex),'Parent',fig2ax2, ...
             'DisplayName','time vs key', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);

set(fig2ax2,'FontSize',14);

ylabel('obs count');
xlabel('time');
datetick(fig2ax2,'x');

%% Create axes for observation index VS. linked list key

fig2ax3 = axes('Parent',figure2,'Position',[0.13 0.11 0.78 0.22]);
box(fig2ax3,'on');

xstring = sprintf('obsmat(:,%d)',obs.indindex);
ystring = sprintf('obsmat(:,%d)',obs.keyindex);
scatter(obsmat(:,obs.indindex),obsmat(:,obs.keyindex),'Parent',fig2ax3, ...
             'DisplayName','count vs key', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);

set(fig2ax3,'FontSize',14);

xlabel('obs count');
ylabel('key');

refreshdata
linkdata on

%% ------------------------------------------------------------------------
% Create figure for 2D scatterplots plots
% -------------------------------------------------------------------------

figure3 = figure(3); clf(figure3); orient tall; wysiwyg


%% Create axes for observation vs ensemble
% This figure is most useful when all the 'bad' obs have been
% replaced by Matlab's NAN so as not to blow the scale.
% The idea is - if both copies 'match', they line up on the diagonal.

fig3ax2 = axes('Parent',figure3,'Position',[0.15 0.05 0.7 0.5]);
xstring = sprintf('obsmat(:,%d)',obs.obsindex);
ystring = sprintf('obsmat(:,%d)',obs.copyindex);

h2 = scatter(obsmat(:,obs.obsindex), obsmat(:,obs.copyindex), ...
    'Parent',fig3ax2, ...
    'DisplayName','copy1 v copy2', ...
    'XDataSource',xstring, ...
    'YDataSource',ystring);

set(fig3ax2,'FontSize',14);

xlabel(obs.colnames{obs.obsindex});
ylabel(obs.colnames{obs.copyindex});
h = title({obs.ObsTypeString, ...
      sprintf('%s ---> %s',obs.timestring(1,:),obs.timestring(2,:)) });
set(h,'Interpreter','none');

axlims = [min(axis) max(axis) min(axis) max(axis)];
axis(axlims)
line([min(axis) max(axis)],[min(axis) max(axis)],'LineWidth',1.5,'Color','k')
grid(fig3ax2,'on');
xmin = axlims(1);
xmax = axlims(2);

if (sum(isfinite(get(h2,'YData')))) == 0
    Print_Empty_Banner(obs.colnames{obs.copyindex});
end

refreshdata
linkdata on

%% Create axes for ObsVal vs. DART QC scatterplot
fig3ax1 = axes('Parent',figure3,'Position',[0.15 0.675 0.7 0.25]);
box(fig3ax1,'on');

xstring = sprintf('obsmat(:,%d)',obs.obsindex);
ystring = sprintf('obsmat(:,%d)',obs.qcindex);

h1 = scatter(obsmat(:,obs.obsindex),obsmat(:,obs.qcindex), ...
    'Parent',fig3ax1, ...
    'DisplayName','obs vs qc', ...
    'XDataSource',xstring, ...
    'YDataSource',ystring);

set(fig3ax1,'FontSize',14);

xlabel(obs.colnames{obs.obsindex});
ylabel(obs.colnames{obs.qcindex});
h = title(obs.ObsTypeString);
set(h,'Interpreter','none');
axis([xmin xmax 0 8]);   % use same horizontal scale as fig3ax2
grid(fig3ax1,'on');


fprintf('QC summary follows:\n')
LabelQC(obs.colnames{obs.qcindex}, obs.qc)

function LabelQC(QCString, qcarray)
%% Create legend for (DART) QC values.
%
% 0     observation assimilated
% 1     observation evaluated only
%   --- everything above this means the prior and posterior are OK
% 2     assimilated, but the posterior forward operator failed
% 3     Evaluated only, but the posterior forward operator failed
%   --- everything above this means only the prior is OK
% 4     prior forward operator failed
% 5     not used
% 6     prior QC rejected
% 7     outlier rejected

dartqc_strings = { 'assimilated', ...
        'observation evaluated only', ...
   'assimilated, post fwd op failed', ...
     'eval only, post fwd op failed', ...
     'prior forward operator failed', ...
                          'not used', ...
                 'prior QC rejected', ...
                  'outlier rejected', ...
                          'reserved'};

switch lower(strtrim(QCString))
   case 'dart quality control',

      qcvals  = unique(qcarray);
      qccount = zeros(size(qcvals));
      for i = 1:length(qcvals)
         qccount(i) = sum(qcarray == qcvals(i));
         fprintf('(%s == %d) %10d obs [%s]\n',strtrim(QCString), qcvals(i), ...
                            qccount(i), dartqc_strings{qcvals(i)+1});
      end

   otherwise,

      qcvals  = unique(qcarray);
      qccount = zeros(size(qcvals));
      for i = 1:length(qcvals)
         qccount(i) = sum(qcarray == qcvals(i));
         fprintf('(%s == %d) %10d obs\n',strtrim(QCString), qcvals(i), qccount(i));
      end
end



function myworldmap(obs)

%%--------------------------------------------------------------------------
% GET THE ELEVATION DATA AND SET UP THE ASSOCIATED COORDINATE DATA
%---------------------------------------------------------------------------

topo = load('topo');     % GET Matlab-native [180x360] ELEVATION DATASET
lats = -89.5:89.5;       % CREATE LAT ARRAY FOR TOPO MATRIX
lons = 0.5:359.5;        % CREATE LON ARRAY FOR TOPO MATRIX
nlon = length(lons);
nlat = length(lats);

%%--------------------------------------------------------------------------
% IF WE NEED TO SWAP HEMISPHERES, DO SO NOW.
% If we didn't explicitly tell it, make a guess.
%---------------------------------------------------------------------------

axis(obs.region);
ax = obs.region;

if (ax(1) < 0)
   lons = lons - 180.0;
   topo.topo = [ topo.topo(:,nlon/2+1:nlon) topo.topo(:,1:nlon/2) ];
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

elev = topo.topo(lat_ind1:lat_ind2,lon_ind1:lon_ind2);
x    = lons(lon_ind1:lon_ind2);
y    = lats(lat_ind1:lat_ind2);

%%--------------------------------------------------------------------------
% Augment the colormap and the CLim so that the lowest color index can be
% forced to a light gray without compromising the data range.

bob      = colormap;
ncolors  = length(bob);
bob(1,:) = 0.7; % set lowest color to be light gray.
colormap(bob);

cmin    = min(obs.obs);
cmax    = max(obs.obs);
dz      = linspace(cmin,cmax,ncolors-1); % desired dynamic range
dclim   = dz(2) - dz(1);
newcmin = cmin - dclim;
clim    = [newcmin cmax]; % add extra bin at bottom that no data will use.
set(gca,'CLim',clim)
colorbar

%%--------------------------------------------------------------------------
% Contour the "subset" - and give the whole thing an appropriate zlevel
% so the continents are either at the top of the plot (for depth), or
% the bottom (for just about everything else.

orgholdstate = ishold;
hold on;

cmd = '[~,h] = contourf(x,y,elev+newcmin,[newcmin newcmin],''k-'');';

switch  lower(obs.Zunits)
    case 'surface'
        set(gca,'Zdir','normal')
        zlevel = ax(5);
    case 'depth'
        set(gca,'Zdir','reverse')
        zlevel = ax(5); % minimum depth
        cmd = '[~,h] = contour(x,y,elev+newcmin,[newcmin newcmin],''k-'');';
    case 'pressure'
        set(gca,'Zdir','reverse')
        zlevel = ax(6); % maximum pressure
    otherwise
        set(gca,'Zdir','normal')
        zlevel = ax(5); % minimum height
end

eval(cmd)

t1    = hgtransform('Parent',gca);
set(h,'Parent',t1);
m     = makehgtform('translate',[0 0 zlevel]);
set(t1,'Matrix',m)

if (orgholdstate == 0), hold off; end;



function y = set_symbol_size(xx)
% seems like a symbol size of 200 is good for a handful of points
% seems like a symbol size of  10 is good for about 5000 points (or more)
% linearly scale between the two ... sucks

if xx > 5000
    y = 10;
elseif xx < 1
    y = 200;
else
    x = [  1 300 500 900  2400 3700 5000];
    z = [200 150 100  60    20   12   10];
    cs = spline(x,[0 z 0]);
    y = round(ppval(cs,xx));
end



function Print_Empty_Banner(mystring)

ax = axis;
dx = mean(ax(1:2));
dy = mean(ax(3:4));
text(dx,dy,sprintf('No "%s" \n with a meaningful DART QC', mystring), ...
    'HorizontalAlignment','center', ...
    'FontSize',20) 

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
