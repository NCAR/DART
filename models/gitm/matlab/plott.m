function p = plott(x,y,xl,yl,p,m,n,s)

%m=50;
%n=round(m/3)
%y should be the same length as x
%x-axis has max resolution of m
%y-axis has max resolution of n
%s is the symbol you want to use for this line
%p is nxm matrix, element (1,1) is displayed in the bottom-left corner

% DART $Id$
% CREDIT: Alexey Morozov

if isnan(xl)
  xl=[min(x) max(x)]
else
  t = (x >= xl(1) & x <= xl(2));
  x=x(t); %discard anything outside the x-scope
  y=y(t);
end

if isnan(yl)
  yl=[min(y) max(y)]
end

d=1:length(x); %each point will be plotted (maybe on top of another)
%d=1:ceil(length(x)/m):length(x); %not sure what it does
xs=round((x(d)-xl(1))/(xl(2)-xl(1))*(m-1)+1);
  ys=round((y(d)-yl(1))/(yl(2)-yl(1))*(n-1)+1);
  s=s*d./d;
s(ys<1 | ys>n)='*'; %highlight the clipped values
ys(ys<1)=1;
ys(ys>n)=n;


%p=sparse(xs,ys,s); %but is sparse
%p = accumarray([xs ys], s); %but adds s if ij is non-unique

ind=sub2ind(size(p),ys,xs);
p(ind)=s; %48+d;

%p=['|'*ones(n,1) flipud(p); '*' '-'*ones(1,m)];
%disp(char(p))

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
