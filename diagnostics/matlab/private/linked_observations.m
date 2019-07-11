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

xlabel(fig1ax1, obs.colnames{obs.lonindex});
ylabel(fig1ax1, obs.colnames{obs.latindex});
zlabel(fig1ax1, obs.colnames{obs.zindex});

tmin = min(obs.time);
tmax = max(obs.time);

pstruct = struct('axis', obs.region, ...
    'clim', [min(obs.obs) max(obs.obs)], ...
    'Ztype', obs.Ztyp(1));
FlatEarth(pstruct);

hb = colorbar;
set(get(hb,'YLabel'),'String',obs.ObsTypeString,'Interpreter','none')

h = title(fig1ax1, {obs.ObsTypeString, ...
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
scatter(obsmat(:,obs.timeindex), obsmat(:,obs.qcindex), ...
    'Parent', fig2ax1, ...
    'DisplayName', 'time vs qc', ...
    'XDataSource', xstring, ...
    'YDataSource', ystring);

set(fig2ax1,'FontSize',14);

ax = axis;
ax(3:4) = [0 8]; %min/max range of DART QC values
axis(ax)
datetick(fig2ax1, 'x')
ylabel(fig2ax1, obs.colnames{obs.qcindex});
grid(fig2ax1,'on');

%% Create axes for observation count/density VS. time

fig2ax2 = axes('Parent',figure2,'Position',[0.13 0.41 0.78 0.22]);
box(fig2ax2,'on');

xstring = sprintf('obsmat(:,%d)',obs.timeindex);
ystring = sprintf('obsmat(:,%d)',obs.indindex);
scatter(obsmat(:,obs.timeindex), obsmat(:,obs.indindex), ...
    'Parent',fig2ax2, ...
    'DisplayName','time vs key', ...
    'XDataSource',xstring, ...
    'YDataSource',ystring);

set(fig2ax2,'FontSize',14);

ylabel(fig2ax2, 'obs count');
xlabel(fig2ax2, 'time');
datetick(fig2ax2, 'x');

%% Create axes for observation index VS. linked list key

fig2ax3 = axes('Parent',figure2,'Position',[0.13 0.11 0.78 0.22]);
box(fig2ax3,'on');

xstring = sprintf('obsmat(:,%d)',obs.indindex);
ystring = sprintf('obsmat(:,%d)',obs.keyindex);
scatter(obsmat(:,obs.indindex), obsmat(:,obs.keyindex), ...
    'Parent',fig2ax3, ...
    'DisplayName','count vs key', ...
    'XDataSource',xstring, ...
    'YDataSource',ystring);

set(fig2ax3,'FontSize',14);

xlabel(fig2ax3, 'obs count');
ylabel(fig2ax3, 'key');

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

xlabel(fig3ax2, obs.colnames{obs.obsindex});
ylabel(fig3ax2, obs.colnames{obs.copyindex});
h = title(fig3ax2, {obs.ObsTypeString, ...
    sprintf('%s ---> %s',obs.timestring(1,:),obs.timestring(2,:)) });
set(h,'Interpreter','none');

% Want to make the graphic 'square'; same values on X,Y axes
axmin = min(axis);
axmax = max(axis);
axis([axmin axmax axmin axmax])
set(fig3ax2,'YTick',get(fig3ax2,'XTick'))

line([axmin axmax],[axmin axmax],'LineWidth',1.5,'Color','k','Parent',fig3ax2)
grid(fig3ax2,'on');

if (sum(isfinite(get(h2,'YData')))) == 0
    Print_Empty_Banner(obs.colnames{obs.copyindex});
end

%% Create axes for ObsVal vs. DART QC scatterplot
fig3ax1 = axes('Parent',figure3,'Position',[0.15 0.675 0.7 0.25]);
box(fig3ax1,'on');

xstring = sprintf('obsmat(:,%d)',obs.obsindex);
ystring = sprintf('obsmat(:,%d)',obs.qcindex);

scatter(obsmat(:,obs.obsindex),obsmat(:,obs.qcindex), ...
    'Parent',fig3ax1, ...
    'DisplayName','obs vs qc', ...
    'XDataSource',xstring, ...
    'YDataSource',ystring);

set(fig3ax1,'FontSize',14);

xlabel(fig3ax1, obs.colnames{obs.obsindex});
ylabel(fig3ax1, obs.colnames{obs.qcindex});
h = title(fig3ax1, obs.ObsTypeString);
set(h,'Interpreter','none');
axis([axmin axmax 0 8]);   % use same horizontal scale as fig3ax2
grid(fig3ax1,'on');

fprintf('\nQC summary follows:\n')
LabelQC(obs.colnames{obs.qcindex}, obs.qc)

refreshdata
linkdata on


%% ------------------------------------------------------------------------


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
    'vertical conversion failed', ...
    'reserved'};

switch lower(strtrim(QCString))
    case 'dart quality control'
        
        qcvals  = unique(qcarray);
        qccount = zeros(size(qcvals));
        for i = 1:length(qcvals)
            qccount(i) = sum(qcarray == qcvals(i));
            fprintf('(%s == %d) %10d obs [%s]\n',strtrim(QCString), qcvals(i), ...
                qccount(i), dartqc_strings{qcvals(i)+1});
        end
        
    otherwise
        
        qcvals  = unique(qcarray);
        qccount = zeros(size(qcvals));
        for i = 1:length(qcvals)
            qccount(i) = sum(qcarray == qcvals(i));
            fprintf('(%s == %d) %10d obs\n',strtrim(QCString), qcvals(i), qccount(i));
        end
end


%% ------------------------------------------------------------------------


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


%% ------------------------------------------------------------------------


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
