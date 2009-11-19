function linked_observations(obsmat,obs)
% linked_observations(obs)
%
% obs is a structure with the following required components
%
% obs.lons	longitudes of the observations
% obs.lats	latitudes of the observations
% obs.z		vertical level (depth) of the observations
% obs.obs	observation values
% obs.qc	observation DART QC code

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Create figure
%figure1 = figure('XVisual',...
%    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
%    'Renderer','OpenGL');
figure1 = figure(1); clf(figure1);;

%% Create axes for 3D plot
axes0 = axes('Parent',figure1,'OuterPosition',[0 0 1 0.90],'FontSize',12);
view(axes0,[-37.5 30]);
grid(axes0,'on');
hold(axes0,'all');

xstring = sprintf('obsmat(:,%d)',obs.lonindex);
ystring = sprintf('obsmat(:,%d)',obs.latindex);
zstring = sprintf('obsmat(:,%d)',obs.zindex  );

h0 = scatter3(obsmat(:,obs.lonindex), obsmat(:,obs.latindex), obsmat(:,obs.zindex), ...
             'Parent',axes0,'DisplayName','observation locations', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring, ...
             'ZDataSource',zstring);

worldmap('light');
xlabel('longitude');
ylabel('latitude');
zlabel('depth');
h = title({obs.ObsTypeString, ...
      sprintf('%s ---> %s',obs.timestring(1,:),obs.timestring(2,:)) });
set(h,'Interpreter','none')
linkdata on

%% Create figure for ancillary plots

figure2 = figure(2); clf(figure2); orient tall; wysiwyg

%% Create axes for time VS. QC
axes4 = axes('Parent',figure2,'OuterPosition',[0 0.80 1 0.175]);
set(axes4,'XAxisLocation','top')
box(axes4,'on');
hold(axes4,'all');

xstring = sprintf('obsmat(:,%d)',obs.timeindex);
ystring = sprintf('obsmat(:,%d)',obs.qcindex);
h4 = scatter(obsmat(:,obs.timeindex),obsmat(:,obs.qcindex),'Parent',axes4, ...
             'DisplayName','time vs qc', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
datetick(axes4,'x',6);
ylabel(obs.QCString);

%% Create axes for observation index VS. time 
axes3 = axes('Parent',figure2,'OuterPosition',[0 0.575 1 0.175]);
box(axes3,'on');
hold(axes3,'all');

xstring = sprintf('obsmat(:,%d)',obs.timeindex);
ystring = sprintf('obsmat(:,%d)',obs.indindex);
h3 = scatter(obsmat(:,obs.timeindex),obsmat(:,obs.indindex),'Parent',axes3, ...
             'DisplayName','time vs key', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
ylabel('obs count');
datetick(axes3,'x',6);

%% Create axes for observation index VS. linked list key
axes2 = axes('Parent',figure2,'OuterPosition',[0.0 0.400 1 0.15]);
box(axes2,'on');
hold(axes2,'all');

xstring = sprintf('obsmat(:,%d)',obs.indindex);
ystring = sprintf('obsmat(:,%d)',obs.keyindex);
h2 = scatter(obsmat(:,obs.indindex),obsmat(:,obs.keyindex),'Parent',axes2, ...
             'DisplayName','count vs key', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
xlabel('obs count');
ylabel('key');

%% Create axes for QC vs. ObsVal scatterplot
axes1 = axes('Parent',figure2,'Position',[0.05 0.05 0.6 0.25]);
box(axes1,'on');
hold(axes1,'all');

xstring = sprintf('obsmat(:,%d)',obs.obsindex);
ystring = sprintf('obsmat(:,%d)',obs.qcindex);
h1 = scatter(obsmat(:,obs.obsindex),obsmat(:,obs.qcindex),'Parent',axes1, ...
             'DisplayName','obs vs qc', ...
             'XDataSource',xstring, ...
             'YDataSource',ystring);
xlabel(obs.CopyString);
title(obs.QCString);


LabelQC(obs.QCString, obs.qc)


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
      for i = 1:length(qcvals)
         qccount(i) = sum(qcarray == qcvals(i));
         s{i} = sprintf('%d - %s %d obs',qcvals(i), dartqc_strings{qcvals(i)+1}, qccount(i));
      end

      set(gca,'YTick',qcvals,'YAxisLocation','right')
      set(gca,'YTickLabel',char(s))

   otherwise,
      str = sprintf('no way to interpret values of %s',strtrim(QCString));
      text(0.0, 0.0, str)
end
