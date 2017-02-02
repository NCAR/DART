%% selerad

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

load sele_obs.jan1-10;
err = sele_obs(:, 1);
N = size(sele_obs,1); % number of locations


%% make a matrix of distances and then find all locations
%  within a certain radius ...

x = sele_obs(:,2);
y = sele_obs(:,3);

X = x * ones(1,N); XT = X';
Y = y * ones(1,N); YT = Y';

d = sqrt( (X -XT).^2 + (Y-YT).^2 );

nmatches = zeros(1,N);
for i=1:N,
   nmatches(i) = sum( d(i,:) < 0.02 );
end

%----------------------------------------------------------------------

counter = 0;
for i = 1:N
  counter = counter + 1;
  xx4(counter) = sele_obs(i, 2)*180.0/3.14159-360;
  yy4(counter) = sele_obs(i, 3)*180.0/3.14159;
end

load coast;

axesm('mercator','MapLatLimit',[10 80], 'MapLonLimit', [-160 -40]); gridm; plotm(lat,long,'color',[.35 .35 .35]); 
tightmap
mh = mlabel;  set(mh,'FontSize',12);   % meridian labels
ph = plabel;  set(ph,'FontSize',12);   % parallel labels

%plotm(yy4, xx4, 'r.', 'markersize', 15);
for i=1:N,
   h = textm(yy4(i),xx4(i),sprintf('%d',nmatches(i)));
   set(h,'HorizontalAlignment','center', ...
         'VerticalAlignment','middle',  ...
         'FontSize',12, ...
         'Fontweight','bold', ...
         'Color',[1 0 0]);
end

%ylabel('Latitude (N)', 'fontsize', 12)
%xlabel('Longitude (E)', 'fontsize', 12)
xlabel('Locations of radiosondes used for verification, June 18-27, 2003', 'fontsize', 14)

 print -dpsc selerad_jun.ps
 print -dpng selerad_jun.png

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
