function h=imagescn(x,y,c)
%IMAGESCN - ImageSC with Nonlinear axis.
%   IMAGESCN(X,Y,C) creates a 2D plot analogous to imagesc in that given
%   the centers of rectangles and colors, it draws them (centered), but
%   DOES allow (as opposed to IMAGESC(X,Y,C)) for nonlinear axis spacing 
%   without throwing away the last row and last column of data, as does
%   PCOLOR(X,Y,C). To do that it has to guess what the rectangle sizes are
%   and particularly it uses LINEAR extrapolation to guess the size of the
%   last row and column.
%
%   Feel free to use as a draft for your own functions.

% DART $Id$
% CREDIT: Alexey Morozov

x=x(:); %collapse any matrices or row vectors into column vectors
y=y(:); %collapse any matrices or row vectors into column vectors

dx=diff(x);
dx=[dx; interp1(1:length(dx),dx,length(dx)+1,'linear','extrap')];
xn=x-dx/2; %the coordinates of new verticies are shifted from centers by the 1/2 rule
xn(end+1)=x(end)+dx(end)/2;

dy=diff(y);
dy=[dy; interp1(1:length(dy),dy,length(dy)+1,'linear','extrap')];
yn=y-dy/2;
yn(end+1)=y(end)+dy(end)/2;

%If you look in "help pcolor", you'll see that pcolor with default shading
%doesn't use the last row and column of c, but needs it to be the same size 
%as x and y, so here, I just copy the values already in there (to keep 
%max(c) and min(c) the same)
c(end+1,:)=c(end,:);
c(:,end+1)=c(:,end);

h=pcolor(xn,yn,c);

% size(xn)
% size(yn)
% size(c)
% clf
% hold on
% for i=1:length(xn)
%     for j=1:length(yn)
%     patch(2*[1 2 2 1]+xn(i), 6*[1 1 2 2]+yn(j),c,'EdgeColor','none')
%     end
% end
% hold off

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
