%% obserr   ... uses output called obserr.dat from somewhere ... 

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%% obsfit.m plot of the error of analysis and guess
clf
a_v=load('obserr.dat');
xa_v=a_v(:,2);
ya_v=a_v(:,1);
xb_v=a_v(:,2);
%
subplot('position', [0.1,0.2,0.7,0.7]), plot(xb_v,ya_v,'k-', xa_v,ya_v,'k--', 'linewidth', 2.0)
axis([0 1.6 0 15])
grid

set(gca, 'FontSize',14)

ylabel('Height (km)', 'fontsize', 14)
xlabel('COSMIC RO excess obs RMS error (%)', 'fontsize', 14);

%
%h = legend('Excess phase', 'Local Ref');
%legpos = get(h,'Position');
%legpos(1:2) = [0.5 0.7];
%set(h,'Position',legpos )
%legend(h, 'boxoff')
%
print -dpsc obserr.ps
print -dpng obserr.png

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
