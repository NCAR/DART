function h2=subview(i,x2,y2,p2)

%i is which plot do you want to subview (1 is the newest added to the
%figure, 2 - the 2nd newest)
%f1 is the figure where to put the subview (assume gcf)
%h1 is the handle to plot to be subvied (includes several "children", like plot itself, legend, colorbar, x,ylabels...)
%h2 is the handle to the subview
%x1 are the xlimits h1 (gotten automatically from h1)
%x2 are the xlimits of the subview (x2 is assumed to be inside x1)
%y2 are ylims of h2 (if set to be nan, then will be picked automatically)
%p2=[.7 .7 .3 .3]; %where and how big to make the subview in h1-relative coordinates

% DART $Id$
% CREDIT: Alexey Morozov

if isnan(p2)
    p2=[.7 .7 .3 .3];
end

f1=gcf;

h1=get(f1,'children');
h1=h1(i);  %assume the plot itself is the last child, discard the rest
% get(h1)
x1=get(h1,'xlim');
y1=get(h1,'ylim');


h2=copyobj(h1,f1);

set(gcf,'CurrentAxes',h2)
set(h2,'FontSize',7)

xlim(x2)

if isnan(y2)
    y2=get(h2,'ylim');
    y2=[max([y1(1) y2(1)]) min([y1(2) y2(2)])];
end
ylim(y2);

xlabel('')
ylabel('')
title('')
p1=get(h1,'position');

p3=[p1(1)+p1(3)*p2(1) p1(2)+p1(4)*p2(2) p1(3)*p2(3) p1(4)*p2(4)]; %-//- in f1-relative coordinates
set(h2,'position', p3);
box on

y2=get(h2,'ylim');

set(gcf,'CurrentAxes',h1)
% r=rectangle('position',[x2(1) y2(1) diff(x2) diff(y2)]);
p4=[p1(1)+p1(3)*abs(x2(1)-x1(1))/diff(x1) p1(2)+p1(4)*abs(y2(1)-y1(1))/diff(y1) p1(3)*diff(x2)/diff(x1) p1(4)*diff(y2)/diff(y1)];
% [p4(1) p4(1)+p4(3) p4(1)+p4(3) p4(1)],[p4(2) p4(2) p4(2)+p4(4) p4(2)+p4(4)]
annotation('line',[p4(1) p4(1)+p4(3)],[p4(2) p4(2)]);
annotation('line',[ p4(1)+p4(3) p4(1)+p4(3) ],[p4(2) p4(2)+p4(4) ]);
annotation('line',[p4(1)+p4(3) p4(1)],[ p4(2)+p4(4) p4(2)+p4(4)]);
annotation('line',[p4(1) p4(1)],[ p4(2)+p4(4) p4(2)]);

if (p4(1)>p3(1) && p4(2)>p3(2)) || (p4(1)<p3(1) && p4(2)<p3(2)) % NE and SW position
    annotation('line',[p4(1) p3(1)],[p4(2)+p4(4) p3(2)+p3(4)]); %line from h1 plot to h2 box
    annotation('line',[p4(1)+p4(3) p3(1)+p3(3)],[p4(2) p3(2)]); %line from h1 plot to h2 box
else % NW and SE position
    annotation('line',[p4(1) p3(1)],[p4(2) p3(2)]); %line from h1 plot to h2 box
    annotation('line',[p4(1)+p4(3) p3(1)+p3(3)],[p4(2)+p4(4) p3(2)+p3(4)]); %line from h1 plot to h2 box
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
