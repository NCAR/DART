
% DART $Id$
% CREDIT: Alexey Morozov

clear
clc
close all
format compact

addpath('~/gitm11/matlab')


cd ../work
ne=3;
nd=14;
aw=1800;

m=112;
n=round(m/3);
p=' '*ones(n,m);

yl=[160 230];
ift=1;

for i=1:ne 
  r=load([ num2str(i) '.dat'])';
  if ift; x=1:length(r); xl=[1 length(r)]; ift=0; end
  p=plott(x,r,xl,yl,p,m,n,i+48);
end

yw=3; %how many characters for ylabel
ya=[num2str(yl(2),['%0' num2str(yw) '.0f']) '|';
' '*ones(n-2,yw) '|'*ones(n-2,1);
num2str(yl(1),['%0' num2str(yw) '.0f']) '|'];

xw=3; %how many characters for ylabel
xa=[' '*ones(1,yw) 'o' '-'*ones(1,m);
' '*ones(1,yw+1) num2str(xl(1),['%0' num2str(xw) '.0f']) ' '*ones(1,m-2*xw) num2str(xl(2),['%0' num2str(xw) '.0f']) ];


p=[ya flipud(p); xa];

disp(' ')
disp(char(p))

cd ../matlab
i=fopen('plot.txt','w');
fwrite(i,[p 10*ones(n+2,1)]');
  fclose(i);

exit

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
