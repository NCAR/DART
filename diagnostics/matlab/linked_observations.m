function linked_observations(obs)
% linked_observations(obs) is a helper function for link_obs.m
% linked_observations is never meant to be called directly.

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

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

% Create figure
%figure1 = figure('XVisual',...
%    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%    'Renderer','OpenGL');
figure1 = figure(1); clf(figure1);

%% Create axes for 3D scatterplot
% should figure out how to query the vertical coordinate to determine
% direction of the Z axis ... 
axes0 = axes('Parent',figure1,'OuterPosition',[0 0 1 0.90],'FontSize',18);
view(axes0,[-37.5 30]);
grid(axes0,'on');
hold(axes0,'all');

xstring = sprintf('obsmat(:,%d)',obs.lonindex);
ystring = sprintf('obsmat(:,%d)',obs.latindex);
zstring = sprintf('obsmat(:,%d)',obs.zindex  );

scatter3(obsmat(:,obs.lonindex), obsmat(:,obs.latindex), obsmat(:,obs.zindex), ...
             'Parent',axes0,'DisplayName','observation locations', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring, ...
             'ZDataSource',zstring);

worldmap('light');

xlabel(obs.colnames{obs.lonindex});
ylabel(obs.colnames{obs.latindex});
zlabel(obs.colnames{obs.zindex});

h = title({obs.ObsTypeString, ...
      sprintf('%s ---> %s',obs.timestring(1,:),obs.timestring(2,:)) });
set(h,'Interpreter','none')
linkdata on

%% Create figure for ancillary plots

figure2 = figure(2); clf(figure2); orient tall; wysiwyg

%% Create axes for time VS. QC
axes4 = axes('Parent',figure2,'OuterPosition',[0 0.80 1 0.175],'FontSize',14);
set(axes4,'XAxisLocation','top')
box(axes4,'on');
hold(axes4,'all');

xstring = sprintf('obsmat(:,%d)',obs.timeindex);
ystring = sprintf('obsmat(:,%d)',obs.qcindex);
scatter(obsmat(:,obs.timeindex),obsmat(:,obs.qcindex),'Parent',axes4, ...
             'DisplayName','time vs qc', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
datetick(axes4,'x',6);
ylabel(obs.colnames{obs.qcindex});

%% Create axes for observation index VS. time 
%axes3 = axes('Parent',figure2,'OuterPosition',[0 0.575 1 0.175]);
axes3 = axes('Parent',figure2,'OuterPosition',[0 0.600 1 0.2],'FontSize',14);
box(axes3,'on');
hold(axes3,'all');

xstring = sprintf('obsmat(:,%d)',obs.timeindex);
ystring = sprintf('obsmat(:,%d)',obs.indindex);
scatter(obsmat(:,obs.timeindex),obsmat(:,obs.indindex),'Parent',axes3, ...
             'DisplayName','time vs key', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
ylabel('obs count');
datetick(axes3,'x',6);

%% Create axes for observation index VS. linked list key
%axes2 = axes('Parent',figure2,'OuterPosition',[0.0 0.375 1 0.2]);
axes2 = axes('Parent',figure2,'OuterPosition',[0.0 0.35 1 0.25],'FontSize',14);
box(axes2,'on');
hold(axes2,'all');

xstring = sprintf('obsmat(:,%d)',obs.indindex);
ystring = sprintf('obsmat(:,%d)',obs.keyindex);
scatter(obsmat(:,obs.indindex),obsmat(:,obs.keyindex),'Parent',axes2, ...
             'DisplayName','count vs key', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
xlabel('obs count');
ylabel('key');

%% Create axes for ObsVal vs. QC scatterplot
axes1 = axes('Parent',figure2,'Position',[0.05 0.05 0.6 0.25],'FontSize',14);
box(axes1,'on');
hold(axes1,'all');

xstring = sprintf('obsmat(:,%d)',obs.obsindex);
ystring = sprintf('obsmat(:,%d)',obs.qcindex);
scatter(obsmat(:,obs.obsindex),obsmat(:,obs.qcindex),'Parent',axes1, ...
             'DisplayName','obs vs qc', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
xlabel( obs.colnames{obs.obsindex});
h = title(  obs.ObsTypeString);
set(h,'Interpreter','none');
LabelQC(obs.colnames{obs.qcindex}, obs.qc)

refreshdata
linkdata on

%% Create axes for observation vs ensemble
% This figure is most useful when all the 'bad' obs have been
% replaced by Matlab's NAN so as not to blow the scale.
% The idea is - if both copies 'match', they line up on the diagonal.

figure3 = figure(3); clf(figure3);

axes5 = axes('Parent',figure3,'OuterPosition',[0 0 1 0.95],'FontSize',18);
grid(axes5,'on');
hold(axes5,'all');

xstring = sprintf('obsmat(:,%d)',obs.obsindex);
ystring = sprintf('obsmat(:,%d)',obs.copyindex);

scatter(obsmat(:,obs.obsindex), obsmat(:,obs.copyindex), ...
             'Parent',axes5,'DisplayName','copy1 v copy2', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);

xlabel(obs.colnames{obs.obsindex});
ylabel(obs.colnames{obs.copyindex});
h = title({obs.ObsTypeString, ...
      sprintf('%s ---> %s',obs.timestring(1,:),obs.timestring(2,:)) });
set(h,'Interpreter','none');

axlims = [min(axis) max(axis) min(axis) max(axis)];
axis(axlims)
plot(axes5,[min(axis) max(axis)],[min(axis) max(axis)],'k-')

%% thats it folks

refreshdata
linkdata on


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
      s = cell(length(qcvals),1);
      for i = 1:length(qcvals)
         qccount(i) = sum(qcarray == qcvals(i));
         s{i} = sprintf('%d - %s %d obs',qcvals(i), dartqc_strings{qcvals(i)+1}, qccount(i));
      end

      set(gca,'YTick',qcvals,'YAxisLocation','right')
      set(gca,'YTickLabel',char(s{:}),'FontSize',12)

   otherwise,
      str = sprintf('no way to interpret values of %s',strtrim(QCString));
      text(0.0, 0.0, str)
end
